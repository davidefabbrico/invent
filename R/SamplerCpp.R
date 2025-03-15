# Function: my_cm
# This function computes the confusion matrix given the predicted labels (sel) and the true labels (true).
# The confusion matrix is a 2x2 matrix that contains the following information:
# - True Positives (TP): Number of cases where the prediction is 1 and the actual value is also 1.
# - True Negatives (TN): Number of cases where the prediction is 0 and the actual value is 0.
# - False Positives (FP): Number of cases where the prediction is 1 but the actual value is 0.
# - False Negatives (FN): Number of cases where the prediction is 0 but the actual value is 1.
#
# Arguments:
# - sel: A vector of predicted labels (either 0 or 1) for each observation.
# - true: A vector of true labels (either 0 or 1) for each observation.
#
# Output:
# - A 2x2 confusion matrix where:
#   - The first row contains [True Negatives, False Positives]
#   - The second row contains [False Negatives, True Positives]

my_cm <- function(sel, true) {
  # Initialize counts for True Positives (tp), True Negatives (tn), False Positives (fp), and False Negatives (fn)
  tp <- 0
  tn <- 0
  fp <- 0
  fn <- 0
  
  # Get the length of the vectors (number of observations)
  p <- length(sel)
  
  # Iterate over each element in the vectors to calculate the confusion matrix components
  for (j in 1:p) {
    # Check for True Positives: predicted and actual are both 1
    if ((sel[j] == true[j]) & (true[j] == 1)) {
      tp <- tp + 1
    }
    # Check for True Negatives: predicted and actual are both 0
    if ((sel[j] == true[j]) & (true[j] == 0)) {
      tn <- tn + 1
    }
    # Check for False Positives: predicted is 1 but actual is 0
    if ((sel[j] != true[j]) & (true[j] == 0)) {
      fp <- fp + 1
    }
    # Check for False Negatives: predicted is 0 but actual is 1
    if ((sel[j] != true[j]) & (true[j] == 1)) {
      fn <- fn + 1
    }
  }
  
  # Create a 2x2 matrix for the confusion matrix
  contTable <- matrix(c(tn, fp, fn, tp), nrow = 2, ncol = 2)
  
  # Return the confusion matrix
  return(contTable)
}

#' invMCMC
#' 
#' @export

##### ------------------------------------------------------------------ ######
invMCMC <- function(y, x, y_val = NULL, x_val = NULL, hyperpar = c(3, 1, 1, 1, 0.00025, 0.4, 1.6, 0.2, 1.8, 0.4, 1.6, 0.2, 1.8), 
                    mht = c(1.4, 0.8, 1, 0.3, 0.7, 0.4, 4, 2.5), 
                    rank = 0.95, iter = 10000, burnin = iter/2, thin = 5, ha = 2, 
                    detail = FALSE, data = NULL, pb = TRUE) {
  
  result <- NULL
  # useful quantities
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
  linearIndex <- c()
  iLin <- 0
  for(j in 1:p) {
    xj <- dx[,j]
    if (all(xj == as.integer(xj))) {
      linearIndex[iLin] <- j
      iLin <- iLin + 1
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
  
  # given that I apply the function lin, the columns are not in the right order
  
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
    # gamma main linear
    gammaLin <- result$gamma_l
    mppi_MainLinear <- apply(gammaLin, 2, mean)
    # gamma main non linear
    gammaNLin <- result$gamma_nl
    mppi_MainNonLinear <- apply(gammaNLin, 2, mean)
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
    mppi_MainLinear <- result$gamma_l
    mppi_MainNonLinear <- result$gamma_nl
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
    true_MainLinear <- as.numeric(data$theta_l[1,] != 0)
    # table linear main effect
    tab_MainLinear <- my_cm(sel_MainLinear, true_MainLinear)
    matt_MainLinear <- mcc(confusionM = tab_MainLinear)
    tpr_MainLinear <- tab_MainLinear[2,2]/sum(tab_MainLinear[,2]) 
    fpr_MainLinear <- tab_MainLinear[2,1]/sum(tab_MainLinear[,1])
    ############### Non linear main effect ###################
    # selected by the model
    sel_MainNonLinear <- as.numeric(mppi_MainNonLinear > 0.5)
    # true non null effect
    true_MainNonLinear <- as.numeric(data$theta_nl[1,] != 0)
    # table linear main effect
    tab_MainNonLinear <- my_cm(sel_MainNonLinear, true_MainNonLinear)
    matt_MainNonLinear <- mcc(confusionM = tab_MainNonLinear)
    if (sum(tab_MainNonLinear[,2]) == 0) {
      tpr_MainNonLinear <- 0
    } else {
      tpr_MainNonLinear <- tab_MainNonLinear[2,2]/sum(tab_MainNonLinear[,2]) 
    }
    if (sum(tab_MainNonLinear[,1]) == 0) {
      fpr_MainNonLinear <- 0
    } else {
      fpr_MainNonLinear <- tab_MainNonLinear[2,1]/sum(tab_MainNonLinear[,1])
    }
    ############### Linear Interaction effect ###################
    sel_IntLinear <- mppi_IntLinear > 0.5
    sel_IntLinear[sel_IntLinear == TRUE] = 1
    true_IntLinear <- data$omega_l != 0
    true_IntLinear[true_IntLinear == TRUE] = 1
    Vec_sel_IntLinear <- sel_IntLinear[upper.tri(sel_IntLinear)]
    Vec_true_IntLinear <- true_IntLinear[upper.tri(true_IntLinear)]
    tab_IntLinear <- my_cm(Vec_sel_IntLinear, Vec_true_IntLinear)
    matt_IntLinear <- mcc(confusionM = tab_IntLinear)
    tpr_IntLinear <- tab_IntLinear[2,2]/sum(tab_IntLinear[,2]) 
    fpr_IntLinear <- tab_IntLinear[2,1]/sum(tab_IntLinear[,1])
    ############### Non linear Interaction effect ###################
    sel_IntNonLinear <- mppi_IntNonLinear > 0.5
    sel_IntNonLinear[sel_IntNonLinear == TRUE] = 1
    TData <- data$omega_nl
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
    tab_IntNonLinear <- my_cm(Vec_sel_IntNonLinear, Vec_true_IntNonLinear)
    matt_IntNonLinear <- mcc(confusionM = tab_IntNonLinear)
    tpr_IntNonLinear <- tab_IntNonLinear[2,2]/sum(tab_IntNonLinear[,2]) 
    fpr_IntNonLinear <- tab_IntNonLinear[2,1]/sum(tab_IntNonLinear[,1])
    # Total metric
    selectInd <- c(as.numeric(sel_MainLinear), as.numeric(sel_MainNonLinear), 
                   Vec_sel_IntLinear , Vec_sel_IntNonLinear)
    trueInd <- c(as.numeric(true_MainLinear), as.numeric(true_MainNonLinear), 
                 Vec_true_IntLinear, Vec_true_IntNonLinear)
    ##########  aggregate
    contTable <- my_cm(selectInd, trueInd)
    matt <- mcc(confusionM = contTable)
    tpr <- contTable[2,2]/sum(contTable[,2]) 
    fpr <- contTable[2,1]/sum(contTable[,1]) 
  }
  # execution time
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
      # if there is no data generating mechanism
      if (!is.null(y_val)) {
        # if there is the posterior predictive distribution
        res <- list(linear_predictor = yhat, y_OutSample = y_tilde, LogLikelihood = ll, y_tildeChain = y_tComplete,
                    mse_inSample = mse_is, mse_outSample = mse_os, mppi_MainLinear = mppi_MainLinear, 
                    mppi_MainNonLinear = mppi_MainNonLinear, mppi_IntLinear = mppi_IntLinear, 
                    mppi_IntNonLinear = mppi_IntNonLinear, Execution_Time = execution_time, acc_rate = acc_rate, sigma = sigma)
      } else {
        res <- list(linear_predictor = yhat, LogLikelihood = ll, 
                    mse = mse_is, mppi_MainLinear = mppi_MainLinear, 
                    mppi_MainNonLinear = mppi_MainNonLinear, mppi_IntLinear = mppi_IntLinear, 
                    mppi_IntNonLinear = mppi_IntNonLinear, Execution_Time = execution_time, acc_rate = acc_rate, sigma = sigma)
      }
    } else {
      if (!is.null(y_val)) {
        # return tpr, fpr, matt
        res <- list(linear_predictor = yhat, y_OutSample = y_tilde, y_tildeChain = y_tComplete,
                    LogLikelihood = ll, mse_inSample = mse_is, 
                    mse_outSample = mse_os, tpr = tpr, tprML = tpr_MainLinear, tprMNL = tpr_MainNonLinear, tprIL = tpr_IntLinear, tprINL = tpr_IntNonLinear, fpr = fpr,
                    fprML = fpr_MainLinear, fprMNL = fpr_MainNonLinear, fprIL = fpr_IntLinear, fprINL = fpr_IntNonLinear,
                    matt = matt, mattML = matt_MainLinear, mattMNL = matt_MainNonLinear, mattIL = matt_IntLinear, mattINL = matt_IntNonLinear,
                    Execution_Time = execution_time, acc_rate = acc_rate, sigma = sigma)
      } else {
        # return tpr, fpr, matt
        res <- list(linear_predictor = yhat, 
                    LogLikelihood = ll, mse = mse_is, tpr = tpr, tprML = tpr_MainLinear, tprMNL = tpr_MainNonLinear, tprIL = tpr_IntLinear, tprINL = tpr_IntNonLinear, fpr = fpr,
                    fprML = fpr_MainLinear, fprMNL = fpr_MainNonLinear, fprIL = fpr_IntLinear, fprINL = fpr_IntNonLinear,
                    matt = matt, mattML = matt_MainLinear, mattMNL = matt_MainNonLinear, mattIL = matt_IntLinear, mattINL = matt_IntNonLinear, Execution_Time = execution_time, acc_rate = acc_rate, sigma = sigma)
      }
    }
  }
  return(res)
}


# Function to compute the Gelman-Rubin R-hat statistic
# This function calculates the R-hat diagnostic to assess 
# the convergence of multiple MCMC chains.
# Input:
# - chains: An array where each column represents an independent 
#   MCMC chain, and each row corresponds to an iteration.
# Output:
# - A single R-hat value. If R-hat is close to 1, it suggests that 
#   the chains have converged to the same distribution.
gelman_rhat <- function(chains) {
  # Number of chains (columns) and iterations (rows)
  m <- dim(chains)[2]  # Number of chains
  n <- dim(chains)[1]  # Number of iterations per chain
  # Compute the mean for each chain
  chain_means <- colMeans(chains)
  # Compute the overall mean across all chains
  global_mean <- mean(chain_means)
  # Compute the between-chain variance (B)
  B <- (n / (m - 1)) * sum((chain_means - global_mean)^2)
  # Compute the within-chain variance (W)
  W <- mean(apply(chains, 2, var))
  # Estimate the marginal variance
  var_hat <- ((n - 1) / n) * W + (1 / n) * B
  # Compute the R-hat statistic
  R_hat <- sqrt(var_hat / W)
  return(R_hat)
}

rhatSinglePar <- function(myres, stringName = "") {
  if (stringName == "") {
    stop("Please provide the name of the parameter to compute the R-hat")
  }
  parameterList <- lapply(myres, `[[`, stringName)
  parameterMatrix <- do.call(cbind, parameterList)
  rhatValue <- gelman_rhat(parameterMatrix)
  return(rhatValue)
}

rhatMainPar <- function(myres, stringName = "") {
  if (stringName == "") {
    stop("Please provide the name of the parameter to compute the R-hat")
  }
  parameterList <- lapply(myres, `[[`, stringName)
  if (stringName == "xi_nl" || stringName == "m_nl") {
    d <- myres[[1]]$d
    q <- sum(d)
    p <- dim(myres[[1]]$X_lin)[2]
    rhatValue <- rep(NA, q)
    for (base in 1:q) {
      parameterMatrix <- sapply(parameterList, function(mat) mat[, base])
      rhatValue[base] <- gelman_rhat(parameterMatrix)
    }
  } else {
    rhatValue <- rep(NA, p)
    for (cov in 1:p) {
      parameterMatrix <- sapply(parameterList, function(mat) mat[, cov])
      rhatValue[cov] <- gelman_rhat(parameterMatrix)
    }
  }
  return(rhatValue)
}

rhatIntPar <- function(myres, stringName = "") {
  if (stringName == "") {
    stop("Please provide the name of the parameter to compute the R-hat")
  }
  p <- dim(myres[[1]]$X_lin)[2]
  d <- myres[[1]]$d
  cd <- c(0, cumsum(d))
  if (stringName == "xi_star_nl" || stringName == "m_star_nl") {
    rhatValue <- c()
    parameterList <- lapply(myres, `[[`, stringName)
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        parameterMatrix <- sapply(parameterList, 
                                  function(lst) sapply(lst, function(mat) mat[j, (cd[k]+1):(cd[k+1])]))
        rhatValue <- c(rhatValue, gelman_rhat(parameterMatrix))
      }
    }
  } else {
    rhatValue <- rep(NA, p*(p-1)/2)
    indComb <- 1
    parameterList <- lapply(myres, `[[`, stringName)
    for (j in 1:(p-1)) {
      for (k in (j+1):p) {
        parameterMatrix <- sapply(parameterList, 
                                     function(lst) sapply(lst, function(mat) mat[j, k]))
        rhatValue[indComb] <- gelman_rhat(parameterMatrix)
        indComb <- indComb + 1
      }
    }
  }
  return(rhatValue)
}

# This function calculates the Effective Sample Size (ESS) of a given sequence of samples. 
# ESS is a measure of the effective number of independent samples in a correlated Markov Chain Monte Carlo (MCMC) chain.
# The function uses autocorrelations for lags up to a specified maximum (`max_lag`) to estimate how much dependence exists
# between the samples and adjusts the total number of samples accordingly. The ESS can be used to assess the efficiency 
# of the MCMC sampling process, with higher ESS values indicating better exploration of the parameter space.
# 
# Parameters:
# - `samples`: A numeric vector of MCMC samples.
# - `max_lag`: The maximum number of lags to consider when calculating autocorrelations (default is 100).
# 
# Returns:
# - A numeric value representing the Effective Sample Size (ESS).
calculate_ESS <- function(samples, max_lag = 100) {
  # Calculate the autocorrelation for lags up to max_lag
  acf_values <- acf(samples, plot = FALSE, lag.max = max_lag)$acf[-1]  # Remove lag 0
  # Sum the autocorrelations up to the maximum lag
  autocorr_sum <- sum(2 * acf_values)
  # Calculate the ESS
  n <- length(samples)
  ESS <- n / (1 + 2 * autocorr_sum)
  return(ESS)
}

#' invParMCMC in Parallel
#' 
#' @export
invParMCMC <- function(y, x, hyperpar = c(3, 1, 1, 1, 0.00025, 0.4, 1.6, 0.2, 1.8, 0.4, 1.6, 0.2, 1.8), 
                      mht = c(1.4, 0.8, 1, 0.3, 0.7, 0.4, 4, 2.5), PredOutOfSample = TRUE,
                      rank = 0.95, iter = 10000, burnin = iter/2, thin = 5, ha = 2, 
                      data = NULL, nchain = 2, percValidationSet = 20, seed = 10) {
  # Check the number of cores
  if (detectCores() < nchain) {
    stop("The number of cores is less than the number of chains")
  } else {
    # Parallel
    nobs <- length(y)
    cor_all <- nchain
    percentage_to_select <- percValidationSet
    num_to_select <- round(nobs * (percentage_to_select / 100))
    registerDoParallel(cores = cor_all)
    registerDoRNG(seed = seed)
    cat("Starting parallel MCMC with", cor_all, "chains, n =", nobs, "and p =", p, "...")
    myres <- foreach(k = 1:cor_all) %dorng% {
      if (PredOutOfSample) {
        random_sample <- sample(1:nobs, num_to_select)
        x_train <- x[-random_sample,]
        y_train <- y[-random_sample]
        x_val <- x[random_sample,]
        y_val <- y[random_sample]
      } else {
        x_val <- NULL
        y_val <- NULL
      }
      # Start the MCMC
      if (is.null(data)) {
        res0 <- invMCMC(y = y_train, x = x_train, y_val = y_val, x_val = x_val, 
                        iter = iter, burnin = burnin, 
                        mht = mht, thin = thin, ha = ha, 
                        hyperpar = hyperpar, detail = TRUE, pb = FALSE)
      } else {
        res0 <- invMCMC(y = y_train, x = x_train, y_val = y_val, x_val = x_val, 
                        iter = iter, burnin = burnin, 
                        mht = mht, thin = thin, ha = ha, 
                        hyperpar = hyperpar, detail = TRUE, pb = FALSE, 
                        data = data)
        # ritornare solo quello che ci interessa
      }
      return(res0)
    } # end parallel
    cat(" Completed!\n\n")
    
    # Compute the R-hat statistic
    # M Star Linear PROBLEMA
    rhatValueMStarLinear <- rhatIntPar(myres, "m_star_l")
    # M Star Non Linear PROBLEMA
    rhatValueMStarNonLinear <- rhatIntPar(myres, "m_star_nl")
    # Theta Linear
    rhatValueThetaLinear <- rhatMainPar(myres, "theta_l")
    # Theta Non Linear
    rhatValueThetaNonLinear <- rhatMainPar(myres, "theta_nl")
    # Xi Linear
    rhatValueXiLinear <- rhatMainPar(myres, "xi_l")
    # Xi Non Linear
    rhatValueXiNonLinear <- rhatMainPar(myres, "xi_nl")
    # M linear
    rhatValueMLinear <- rhatMainPar(myres, "m_l")
    # M Non Linear
    rhatValueMNonLinear <- rhatMainPar(myres, "m_nl")
    # Pi Linear
    rhatValuePiLinear <- rhatSinglePar(myres, "pi_l")
    # Pi Non Linear
    rhatValuePiNonLinear <- rhatSinglePar(myres, "pi_nl")
    # Pi Star Linear
    rhatValuePiStarLinear <- rhatSinglePar(myres, "pi_star_l")
    # Pi Star Non Linear
    rhatValuePiStarNonLinear <- rhatSinglePar(myres, "pi_star_nl")
    # Alpha Star Linear
    rhatValueAlphaStarLinear <- rhatIntPar(myres, "alpha_star_l")
    # Alpha Star Non Linear
    rhatValueAlphaStarNonLinear <- rhatIntPar(myres, "alpha_star_nl")
    # Tau Linear
    rhatValueTauLinear <- rhatMainPar(myres, "tau_l")
    # Tau Non Linear
    rhatValueTauNonLinear <- rhatMainPar(myres, "tau_nl")
    # Tau Star Linear
    rhatValueTauStarLinear <- rhatIntPar(myres, "tau_star_l")
    # Tau Star Non Linear
    rhatValueTauStarNonLinear <- rhatIntPar(myres, "tau_star_nl")
    # Xi Star Linear
    rhatValueXiStarLinear <- rhatIntPar(myres, "xi_star_l")
    # Xi Star Non Linear
    rhatValueXiStarNonLinear <- rhatIntPar(myres, "xi_star_nl")
    # Intercept
    rhatValueIntercept <- rhatSinglePar(myres, "intercept")
    # Model Variance
    rhatValueSigma <- rhatSinglePar(myres, "sigma")
    
    # Collect all Rhat values into a vector
    all_rhat <- c(
      # Interaction parameters (matrix/array parameters)
      rhatValueMStarLinear, rhatValueMStarNonLinear,
      
      # Main effect parameters (vector parameters)
      rhatValueThetaLinear, rhatValueThetaNonLinear,
      rhatValueXiLinear, rhatValueXiNonLinear,
      rhatValueMLinear, rhatValueMNonLinear,
      
      # Probability parameters (scalars)
      rhatValuePiLinear, rhatValuePiNonLinear,
      rhatValuePiStarLinear, rhatValuePiStarNonLinear,
      
      # Interaction coefficients
      rhatValueAlphaStarLinear, rhatValueAlphaStarNonLinear,
      
      # Variance parameters
      rhatValueTauLinear, rhatValueTauNonLinear,
      rhatValueTauStarLinear, rhatValueTauStarNonLinear,
      
      # Interaction weights
      rhatValueXiStarLinear, rhatValueXiStarNonLinear,
      
      # Model fundamentals
      rhatValueIntercept, rhatValueSigma
    )
    
    # 1. Check for NA/NaN values 
    if (anyNA(all_rhat) || any(is.nan(all_rhat))) {
      warning("WARNING: NA/NaN values detected in Rhat statistics!")
      message("Problematic positions:")
      print(which(is.na(all_rhat) | is.nan(all_rhat)))
    }
    
    # 2. Convergence check (Rhat <= 1.1 threshold)
    threshold <- 1.1  # Standard convergence threshold
    valid_rhat <- all_rhat[!is.na(all_rhat) & !is.nan(all_rhat)]
    good_rhat_count <- sum(valid_rhat <= threshold, na.rm = TRUE)
    total_valid <- length(valid_rhat)
    
    if (total_valid == 0) {
      stop("ERROR: No valid Rhat values available for analysis")
    }
    
    convergence_ratio <- good_rhat_count / total_valid
    
    if (convergence_ratio >= 0.95) {
      message(sprintf(
        "\u2705 %.1f%% of Rhat values (%d/%d) <= %.2f",
        convergence_ratio * 100,
        good_rhat_count,
        total_valid,
        threshold
      ))
    } else {
      warning(sprintf(
        "\u274c WARNING: Only %.1f%% of Rhat values (n = %d/%d) <= %.2f",
        convergence_ratio * 100,
        good_rhat_count,
        total_valid,
        threshold
      ))
      
      # Show problematic values
      message("Non-converged parameters (Rhat > ", threshold, "):")
      print(round(valid_rhat[valid_rhat > threshold], 3))
    }
    
    # 3. Detailed diagnostics report 
    message("\nRhat Statistics Summary:")
    print(summary(valid_rhat))
    
    # Compute the ESS
    
  } # end check cores
  
  return(myres)
}


# my_model <- function(formula, data) {
#   X <- model.matrix(formula, data)
#   y <- model.response(model.frame(formula, data))
# }






