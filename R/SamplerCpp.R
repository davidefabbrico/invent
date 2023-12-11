#' invMCMC
#' 
#' @export

##### ------------------------------------------------------------------ ######
invMCMC <- function(y, x, y_val = NULL, x_val = NULL, hyperpar = c(5, 2, 5, 5, 0.00025, 0.4, 1.6, 0.2, 1.8, 0.4, 1.6, 0.2, 1.8), 
                   mht = c(1.4, 0.8, 1, 0.3, 0.7, 0.4, 4, 2.5), 
                   rank = 0.95, iter = 10000, burnin = iter/2, thin = 5, ha = 2, 
                   detail = FALSE, data = NULL, pb = TRUE) {
  
  result <- NULL
  #useful quantities
  nobs <- dim(x)[1] # number of observations
  p <- dim(x)[2] # number of covariates
  n_cat <- 0
  nval <- 0
  if (!is.null(x_val)) {
    nval <- dim(x_val)[1]
  }
  # matrice del disegno completa sia training che validation
  dx <- rbind(x, x_val)
  #build our basis p spline representation
  X <- vector()
  X_lin <- vector()
  X_nlin <- vector()
  d <- vector()
  
  inD <- 1
  for(j in 1:p) {
    xj <- dx[,j]
    if (any(xj != as.integer(xj))) {
      X <- cbind(X, xj)
      xjl <- lin(xj)
      X_lin <- cbind(X_lin, xjl)
      xjtilde <- sm(x = xj, rankZ = rank) 
      X_nlin <- cbind(X_nlin, xjtilde)
      d[inD] <- dim(xjtilde)[2]
      inD <- inD + 1
    }       
  }
  # linear covariates transformation
  for(j in 1:p) {
    xj <- dx[,j]
    if (all(xj == as.integer(xj))) {
      X <- cbind(X, xj)
      n_cat <- n_cat + 1
      X_lin <- cbind(X_lin, xj)
    }   
  }
  cd <- c(0, cumsum(d))
  q <- sum(d)
  
  if (nval != 0) {
    X_l <- X_lin[1:nobs,]
    X_val_l <- X_lin[(nobs+1):(nobs+nval),]
    X_nl <- X_nlin[1:nobs,]
    X_val_nl <- X_nlin[(nobs+1):(nobs+nval),]
  } else {
    X_l <- X_lin
    X_val_l <- matrix(0, nrow = 0, ncol = 0)
    X_nl <- X_nlin
    X_val_nl <- matrix(0, nrow = 0, ncol = 0)
  }
  
  # Prediction
  # build our basis p spline representation
  # X_val <- vector()
  # X_val_l <- vector()
  # X_val_nl <- vector()
  # d_val <- vector()
  # 
  # if (!is.null(y_val)) {
  #   inD_val <- 1
  #   for(j in 1:p) {
  #     xj_val <- x_val[,j]
  #     if (any(xj_val != as.integer(xj_val))) {
  #       X_val <- cbind(X_val, xj_val)
  #       xjl_val <- lin(xj_val)
  #       X_val_l <- cbind(X_val_l, xjl_val)
  #       xjtilde_val <- sm(x = xj_val, rankZ = rank) 
  #       X_val_nl <- cbind(X_val_nl, xjtilde_val)
  #       d_val[inD_val] <- dim(xjtilde_val)[2]
  #       inD_val <- inD_val + 1
  #     }       
  #   }
  #   inD <- 1
  #   for(j in 1:p) {
  #     xj <- x[,j]
  #     if (any(xj != as.integer(xj))) {
  #       X <- cbind(X, xj)
  #       xjl <- lin(xj)
  #       X_l <- cbind(X_l, xjl)
  #       xjtilde <- sm(x = xj, rankZ = rank) 
  #       X_nl <- cbind(X_nl, xjtilde)
  #       d[inD] <- dim(xjtilde)[2]
  #       inD <- inD + 1
  #     }       
  #   }
  #   # linear covariates transformation
  #   for(j in 1:p) {
  #     xj_val <- x_val[,j]
  #     if (all(xj_val == as.integer(xj_val))) {
  #       X_val <- cbind(X_val, xj_val)
  #       # n_cat <- n_cat + 1
  #       X_val_l <- cbind(X_val_l, xj_val)
  #     }   
  #   }
  #   cd_val <- c(0, cumsum(d_val))
  #   q_val <- sum(d_val)
  # } else {
  #   cd_val <- 0
  #   q_val <- 0
  # }
  
  # Call the C++ function
  result = bodyMCMC(as.vector(y), as.integer(p), as.integer(nobs), as.vector(cd),
                    as.vector(d), as.matrix(X_l), as.matrix(X_nl), as.matrix(X_val_l), as.matrix(X_val_nl),
                    as.vector(hyperpar), as.vector(mht), as.integer(n_cat), as.integer(iter), 
                    as.integer(burnin), as.integer(thin), ha, as.logical(detail), as.logical(pb))
  
  nlp <- p - n_cat
  nout <- (iter-burnin)/thin
  if (detail == TRUE) {
    # Compute the main metrics
    gammaStarLin <- array(unlist(result$gamma_star_l), dim = c(p, p, nout))
    gammaStarNLin <- array(unlist(result$gamma_star_nl), dim = c(nlp, nlp, nout))
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
  } else {
    mppi_IntLinear <- result$gamma_star_l
    mppi_IntNonLinear <- result$gamma_star_nl
    mppi_MainLinear <- result$gamma_0_l
    mppi_MainNonLinear <- result$gamma_0_nl
    yhat <- result$linear_predictor
  }
  # mean square error
  mse_is <- mse(yhat, y)
  # log-likelihood
  ll <- result$LogLikelihood
  # prediction
  if (!is.null(y_val)) {
    if (detail == TRUE) {
      lp_os <- result$y_oos
      y_tilde <- apply(lp_os, 2, mean) 
    } else {
      y_tilde <- result$y_oos
    }
    mse_os <- mse(y_tilde, y_val)
  } else {
    lp_os <- NULL
    y_tilde <- NULL
    mse_os <- NULL
  }
  if (!is.null(data)) {
    ############### Linear main effect ###################
    # selected by the model
    sel_MainLinear <- as.numeric(mppi_MainLinear > 0.5)
    # true non null effect
    true_MainLinear <- as.numeric(data$alpha_0_l[1,] != 0)
    # table linear main effect
    tab_MainLinear <- matrix(table(sel_MainLinear, true_MainLinear), nrow = 2, ncol = 2)
    matt_MainLinear <- mcc(confusionM = tab_MainLinear)
    tpr_MainLinear <- tab_MainLinear[2,2]/sum(tab_MainLinear[,2]) 
    fpr_MainLinear <- tab_MainLinear[2,1]/sum(tab_MainLinear[,1])
    ############### Non linear main effect ###################
    # selected by the model
    sel_MainNonLinear <- as.numeric(mppi_MainNonLinear > 0.5)
    # true non null effect
    true_MainNonLinear <- as.numeric(data$alpha_0_tilde[1,] != 0)
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
    selectInd <- c(as.numeric(sel_MainLinear), as.numeric(sel_MainNonLinear), 
                   Vec_sel_IntLinear , Vec_sel_IntNonLinear)
    trueInd <- c(as.numeric(true_MainLinear), as.numeric(true_MainNonLinear), 
                 Vec_true_IntLinear, Vec_true_IntNonLinear)
    ##########  aggregate
    contTable <- matrix(table(selectInd, trueInd), nrow = 2, ncol = 2)
    matt <- mcc(confusionM = contTable)
    tpr <- contTable[2,2]/sum(contTable[,2]) 
    fpr <- contTable[2,1]/sum(contTable[,1]) 
  }
  # exexution time
  execution_time <- result$Execution_Time
  # acceptance ratio 
  acc_rate <- result$acc_rate
  # sigma error variance
  sigma <- result$sigma
  # complete y_tilde
  y_tComplete <- result$y_tilde
  # Choice to have the details or not
  if (detail == TRUE) {
    res <- result
  } else {
    if (is.null(data)) {
      if (!is.null(y_val)) {
        res <- list(linear_predictor = yhat, y_OutSample = y_tilde, LogLikelihood = ll, y_tComp = y_tComplete,
                    mse_inSample = mse_is, mse_outSample = mse_os, mppi_MainLinear = mppi_MainLinear, 
                    mppi_MainNonLinear = mppi_MainNonLinear, mppi_IntLinear = mppi_IntLinear, 
                    mppi_IntNonLinear = mppi_IntNonLinear, execution_time = execution_time, acc_rate = acc_rate, sigma = sigma)
      } else {
        res <- list(linear_predictor = yhat, LogLikelihood = ll, 
                    mse = mse_is, mppi_MainLinear = mppi_MainLinear, 
                    mppi_MainNonLinear = mppi_MainNonLinear, mppi_IntLinear = mppi_IntLinear, 
                    mppi_IntNonLinear = mppi_IntNonLinear, execution_time = execution_time, acc_rate = acc_rate, sigma = sigma)
      }
    } else {
      if (!is.null(y_val)) {
        # return tpr, fpr, matt
        res <- list(linear_predictor = yhat, y_OutSample = y_tilde, y_tComp = y_tComplete,
                    LogLikelihood = ll, mse_inSample = mse_is, 
                    mse_outSample = mse_os, tpr = tpr, fpr = fpr,
                    matt = matt, execution_time = execution_time, acc_rate = acc_rate, sigma = sigma)
      } else {
        # return tpr, fpr, matt
        res <- list(linear_predictor = yhat, 
                    LogLikelihood = ll, mse = mse_is, tpr = tpr, fpr = fpr,
                    matt = matt, execution_time = execution_time, acc_rate = acc_rate, sigma = sigma)
      }
    }
  }
  
  return(res)
}
