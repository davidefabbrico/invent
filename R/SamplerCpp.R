#' invMCMC
#' 
#' @export

##### ------------------------------------------------------------------ ######
invMCMC <- function(y, x, hyperpar = c(5, 25, 5, 5, 0.00025, 0.4, 1.6, 0.2, 1.8), 
                   mht = c(1.4, 0.8, 1, 0.3, 0.7, 0.4, 4, 2.5), 
                   rank = 0.95, iter = 10000, burnin = iter/2, thin = 5, ha = 2, n_val = 100, pred = TRUE, 
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
  if (pred) {
    # y
    y_val <- y[(nobs-n_val+1):nobs]
    y <- y[1:(nobs-n_val)]
    # x
    x_val <- cbind(X_l[(nobs-n_val+1):nobs,], X_nl[(nobs-n_val+1):nobs,])
    X_l <- X_l[1:(nobs-n_val),]
    X_nl <- X_nl[1:(nobs-n_val),]
    nobs = nobs - n_val
  } else {
    x_val <- matrix(0, nrow = n_val, ncol = p+q)
  }
  result = bodyMCMC(as.vector(y), as.integer(p), as.integer(nobs), as.vector(cd), 
                    as.vector(d), as.matrix(X_l), as.matrix(X_nl), 
                    as.vector(hyperpar), as.vector(mht), as.integer(iter), 
                    as.integer(burnin), as.integer(thin), ha, as.matrix(x_val), as.integer(pred))
  
  if (is.null(data)) {
    print("End of MCMC... MPPI Evaluation")
  } else {
    print("End of MCMC... Metric Evaluation")
  }
  nout <- (iter-burnin)/thin
  # Compute the main metrics
  gammaStarLin <- array(unlist(result$gamma_star_l), dim = c(p, p, nout))
  gammaStarNLin <- array(unlist(result$gamma_star_nl), dim = c(p, p, nout))
  # gamma 0 linear
  gamma0Lin <- result$gamma_0_l[-c(1:burnin),]
  mppi_MainLinear <- apply(gamma0Lin, 2, mean)
  # gamma 0 non linear
  gamma0NLin <- result$gamma_0_nl[-c(1:burnin),]
  mppi_MainNonLinear <- apply(gamma0NLin, 2, mean)
  # gamma star linear (list of matrix)
  gammaStarLin <- gammaStarLin[,,-c(1:burnin)]
  mppi_IntLinear <- apply(gammaStarLin, c(1,2), mean)
  # gamma star non linear
  gammaStarNLin <- gammaStarNLin[,,-c(1:burnin)]
  mppi_IntNonLinear <- apply(gammaStarNLin, c(1,2), mean)
  # linear predictor
  lp_is <- result$linear_predictor[-c(1:burnin),]
  yhat <- apply(lp_is, 2, mean)
  # mean square error
  mse_is <- mse(yhat, y)
  # log-likelihood
  ll <- result$LogLikelihood
  # prediction
  if (pred) {
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
    tpr_MainLinear <- tab_MainLinear[1,1]/sum(tab_MainLinear[,1]) 
    fpr_MainLinear <- tab_MainLinear[1,2]/sum(tab_MainLinear[,2])
    ############### Non linear main effect ###################
    # selected by the model
    sel_MainNonLinear <- factor(as.numeric(mppi_MainNonLinear > 0.5), levels = c(0,1))
    # true non null effect
    true_MainNonLinear <- factor(as.numeric(data$alpha_0_tilde[1,] != 0), levels = c(0,1))
    # table linear main effect
    tab_MainNonLinear <- matrix(table(sel_MainNonLinear, true_MainNonLinear), nrow = 2, ncol = 2)
    matt_MainNonLinear <- mcc(confusionM = tab_MainNonLinear)
    tpr_MainNonLinear <- tab_MainNonLinear[1,1]/sum(tab_MainNonLinear[,1]) 
    fpr_MainNonLinear <- tab_MainNonLinear[1,2]/sum(tab_MainNonLinear[,2])
    ############### Linear Interaction effect ###################
    sel_IntLinear <- mppi_IntLinear > 0.5
    sel_IntLinear[sel_IntLinear == TRUE] = 1
    true_IntLinear <- data$omega_l != 0
    true_IntLinear[true_IntLinear == TRUE] = 1
    Vec_sel_IntLinear <- sel_IntLinear[upper.tri(sel_IntLinear)]
    Vec_true_IntLinear <- true_IntLinear[upper.tri(true_IntLinear)]
    tab_IntLinear <- matrix(table(Vec_sel_IntLinear, Vec_true_IntLinear), nrow = 2, ncol = 2)
    matt_IntLinear <- mcc(confusionM = tab_IntLinear)
    tpr_IntLinear <- tab_IntLinear[1,1]/sum(tab_IntLinear[,1]) 
    fpr_IntLinear <- tab_IntLinear[1,2]/sum(tab_IntLinear[,2])
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
    tpr_IntNonLinear <- tab_IntNonLinear[1,1]/sum(tab_IntNonLinear[,1]) 
    fpr_IntNonLinear <- tab_IntNonLinear[1,2]/sum(tab_IntNonLinear[,2])
    # Total metric
    selectInd <- c(as.numeric(sel_MainLinear)-1, as.numeric(sel_MainNonLinear)-1, Vec_sel_IntLinear , Vec_sel_IntNonLinear)
    trueInd <- c(as.numeric(true_MainLinear)-1, as.numeric(true_MainNonLinear)-1, Vec_true_IntLinear, Vec_true_IntNonLinear)
    ##########  aggregate
    contTable <- matrix(table(selectInd, trueInd), nrow = 2, ncol = 2)
    matt <- mcc(confusionM = contTable)
    tpr <- contTable[1,1]/sum(contTable[,1]) 
    fpr <- contTable[1,2]/sum(contTable[,2]) 
  }
  # exexution time
  execution_time <- result$Execution_Time
  
  # Choice to have the details or not
  if (detail == TRUE) {
    res <- result
  } else {
    if (is.null(data)) {
      res <- list(linear_predictor = yhat, y_OutSample = y_tilde, LogLikelihood = ll, mse_inSample = mse_is, 
                  mse_outSample = mse_os, mppi_MainLinear = mppi_MainLinear, mppi_MainNonLinear = mppi_MainNonLinear, 
                  mppi_IntLinear = mppi_IntLinear, mppi_IntNonLinear = mppi_IntNonLinear, 
                  execution_time = execution_time)
    } else {
      # return tpr, fpr, matt
      res <- list(linear_predictor = yhat, y_OutSample = y_tilde, 
                  LogLikelihood = ll, mse_inSample = mse_is, 
                  mse_outSample = mse_os, tpr = tpr, fpr = fpr,
                  matt = matt, execution_time = execution_time)
    }
  }
  
  return(res)
}
