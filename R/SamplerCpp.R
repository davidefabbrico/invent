#' invMCMC
#' 
#' @export

##### ------------------------------------------------------------------ ######
invMCMC <- function(y, x, y_val = NULL, x_val = NULL, hyperpar = c(5, 25, 5, 5, 0.00025, 0.4, 1.6, 0.2, 1.8, 0.4, 1.6, 0.2, 1.8), 
                   mht = c(1.4, 0.8, 1, 0.3, 0.7, 0.4, 4, 2.5), 
                   rank = 0.95, iter = 10000, burnin = iter/2, thin = 5, ha = 2, 
                   detail = FALSE, data = NULL) {
  
  result <- NULL
  #useful quantities
  nobs <- dim(x)[1] # number of observations
  p <- dim(x)[2] # number of covariates
  
  #build our basis p spline representation
  X <- vector()
  X_l <- vector()
  X_nl <- vector()
  d <- vector()
  
  for(j in 1:p) {
    xj <- x[,j]
    X <- cbind(X, xj)
    xjl <- lin(xj)
    X_l <- cbind(X_l, xjl)
    xjtilde <- sm(x = xj, rankZ = rank) 
    X_nl <- cbind(X_nl, xjtilde)
    d[j] <- dim(xjtilde)[2]
  }
  cd <- c(0, cumsum(d))
  q <- sum(d)

  # Prediction
  # build our basis p spline representation
  X_val <- vector()
  X_val_l <- vector()
  X_val_nl <- vector()
  d_val <- vector()
  cd_val <- 0
  q_val <- 0
  if (!is.null(y_val)) {
    for(j in 1:p) {
      xj_val <- x_val[,j]
      X_val <- cbind(X_val, xj_val)
      xjl_val <- lin(xj_val)
      X_val_l <- cbind(X_val_l, xjl_val)
      xjtilde_val <- sm(x = xj_val, rankZ = rank) 
      X_val_nl <- cbind(X_val_nl, xjtilde_val)
      d_val[j] <- dim(xjtilde_val)[2]
    }
    cd_val <- c(0, cumsum(d_val))
    q_val <- sum(d_val)
  }
  
  result = bodyMCMC(as.vector(y), as.integer(p), as.integer(nobs), as.vector(cd), as.vector(cd_val),
                    as.vector(d), as.vector(d_val), as.matrix(X_l), as.matrix(X_nl), as.matrix(X_val_l), as.matrix(X_val_nl),
                    as.vector(hyperpar), as.vector(mht), as.integer(iter), 
                    as.integer(burnin), as.integer(thin), ha)
  
  nout <- (iter-burnin)/thin
  # Compute the main metrics
  gammaStarLin <- array(unlist(result$gamma_star_l), dim = c(p, p, nout))
  gammaStarNLin <- array(unlist(result$gamma_star_nl), dim = c(p, p, nout))
  # gamma 0 linear
  gamma0Lin <- result$gamma_0_l
  mppi_MainLinear <- apply(gamma0Lin, 2, mean)
  # gamma 0 non linear
  gamma0NLin <- result$gamma_0_nl
  mppi_MainNonLinear <- apply(gamma0NLin, 2, mean)
  # gamma star linear (list of matrix)
  mppi_IntLinear <- apply(gammaStarLin, c(1,2), mean)
  # gamma star non linear
  mppi_IntNonLinear <- apply(gammaStarNLin, c(1,2), mean)
  # linear predictor
  lp_is <- result$linear_predictor
  yhat <- apply(lp_is, 2, mean)
  # mean square error
  mse_is <- mse(yhat, y)
  # log-likelihood
  ll <- result$LogLikelihood
  # prediction
  if (!is.null(y_val)) {
    lp_os <- result$y_oos
    y_tilde <- apply(lp_os, 2, mean)
    mse_os <- mse(y_tilde, y_val)
  } else {
    lp_os <- NULL
    y_tilde <- NULL
    mse_os <- NULL
  }
  if (!is.null(data)) {
    ############### Linear main effect ###################
    # selected by the model
    sel_MainLinear <- factor(as.numeric(mppi_MainLinear > 0.5), levels = c(0,1))
    # true non null effect
    true_MainLinear <- factor(as.numeric(data$alpha_0_l[1,] != 0), levels = c(0,1))
    # table linear main effect
    tab_MainLinear <- matrix(table(sel_MainLinear, true_MainLinear), nrow = 2, ncol = 2)
    matt_MainLinear <- mcc(confusionM = tab_MainLinear)
    tpr_MainLinear <- tab_MainLinear[2,2]/sum(tab_MainLinear[,2]) 
    fpr_MainLinear <- tab_MainLinear[2,1]/sum(tab_MainLinear[,1])
    ############### Non linear main effect ###################
    # selected by the model
    sel_MainNonLinear <- factor(as.numeric(mppi_MainNonLinear > 0.5), levels = c(0,1))
    # true non null effect
    true_MainNonLinear <- factor(as.numeric(data$alpha_0_tilde[1,] != 0), levels = c(0,1))
    # table linear main effect
    tab_MainNonLinear <- matrix(table(sel_MainNonLinear, true_MainNonLinear), nrow = 2, ncol = 2)
    matt_MainNonLinear <- mcc(confusionM = tab_MainNonLinear)
    tpr_MainNonLinear <- tab_MainNonLinear[2,2]/sum(tab_MainNonLinear[,2]) 
    fpr_MainNonLinear <- tab_MainNonLinear[2,1]/sum(tab_MainNonLinear[,1])
    ############### Linear Interaction effect ###################
    sel_IntLinear <- mppi_IntLinear > 0.5
    sel_IntLinear[sel_IntLinear == TRUE] = 1
    true_IntLinear <- data$omega_l != 0
    true_IntLinear[true_IntLinear == TRUE] = 1
    Vec_sel_IntLinear <- sel_IntLinear[upper.tri(sel_IntLinear)]
    Vec_true_IntLinear <- true_IntLinear[upper.tri(true_IntLinear)]
    tab_IntLinear <- matrix(table(Vec_sel_IntLinear, Vec_true_IntLinear), nrow = 2, ncol = 2)
    matt_IntLinear <- mcc(confusionM = tab_IntLinear)
    tpr_IntLinear <- tab_IntLinear[2,2]/sum(tab_IntLinear[,2]) 
    fpr_IntLinear <- tab_IntLinear[2,1]/sum(tab_IntLinear[,1])
    ############### Non linear Interaction effect ###################
    sel_IntNonLinear <- mppi_IntNonLinear > 0.5
    sel_IntNonLinear[sel_IntNonLinear == TRUE] = 1
    TData <- data$omega_tilde
    true_IntNonLinear <- matrix(0, nrow = p, ncol = p)
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        if (sum(TData[i, (cd[j]+1):(cd[j+1])]) != 0) {
          true_IntNonLinear[i,j] <- 1
        }
      } 
    }
    Vec_sel_IntNonLinear <- sel_IntNonLinear[upper.tri(sel_IntNonLinear)]
    Vec_true_IntNonLinear <- true_IntNonLinear[upper.tri(true_IntNonLinear)]
    tab_IntNonLinear <- matrix(table(Vec_sel_IntNonLinear, Vec_true_IntNonLinear), nrow = 2, ncol = 2)
    matt_IntNonLinear <- mcc(confusionM = tab_IntNonLinear)
    tpr_IntNonLinear <- tab_IntNonLinear[2,2]/sum(tab_IntNonLinear[,2]) 
    fpr_IntNonLinear <- tab_IntNonLinear[2,1]/sum(tab_IntNonLinear[,1])
    # Total metric
    selectInd <- c(as.numeric(sel_MainLinear)-1, as.numeric(sel_MainNonLinear)-1, 
                   Vec_sel_IntLinear , Vec_sel_IntNonLinear)
    trueInd <- c(as.numeric(true_MainLinear)-1, as.numeric(true_MainNonLinear)-1, 
                 Vec_true_IntLinear, Vec_true_IntNonLinear)
    ##########  aggregate
    contTable <- matrix(table(selectInd, trueInd), nrow = 2, ncol = 2)
    matt <- mcc(confusionM = contTable)
    tpr <- contTable[2,2]/sum(contTable[,2]) 
    fpr <- contTable[2,1]/sum(contTable[,1]) 
  }
  # exexution time
  execution_time <- result$Execution_Time
  
  # Choice to have the details or not
  if (detail == TRUE) {
    res <- result
  } else {
    if (is.null(data)) {
      if (is.null(y_val)) {
        res <- list(linear_predictor = yhat, y_OutSample = y_tilde, LogLikelihood = ll, 
                    mse_inSample = mse_is, mse_outSample = mse_os, mppi_MainLinear = mppi_MainLinear, 
                    mppi_MainNonLinear = mppi_MainNonLinear, mppi_IntLinear = mppi_IntLinear, 
                    mppi_IntNonLinear = mppi_IntNonLinear, execution_time = execution_time)
      } else {
        res <- list(linear_predictor = yhat, LogLikelihood = ll, 
                    mse = mse_is, mppi_MainLinear = mppi_MainLinear, 
                    mppi_MainNonLinear = mppi_MainNonLinear, mppi_IntLinear = mppi_IntLinear, 
                    mppi_IntNonLinear = mppi_IntNonLinear, execution_time = execution_time)
      }
    } else {
      if (is.null(y_val)) {
        # return tpr, fpr, matt
        res <- list(linear_predictor = yhat, y_OutSample = y_tilde, 
                    LogLikelihood = ll, mse_inSample = mse_is, 
                    mse_outSample = mse_os, tpr = tpr, fpr = fpr,
                    matt = matt, execution_time = execution_time)
      } else {
        # return tpr, fpr, matt
        res <- list(linear_predictor = yhat, 
                    LogLikelihood = ll, mse = mse_is, tpr = tpr, fpr = fpr,
                    matt = matt, execution_time = execution_time)
      }
    }
  }
  
  return(res)
}
