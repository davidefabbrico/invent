##### ------------------------------------------------------------------ ######
# MCMC code implementing non-linear interaction model as varying coefficient  #
# Authors: D. Fabbrico, M. Pedone                                             #  
# utils: file containing useful functions to be called in the main            # 
# function or in the scripts

##### ------------------------------------------------------------------ ######
# Notation
# n_obs: number of observations
# p: number of covariates
#
# d: vector of number of basis used to represent each covariate (p x 1)
# epsilon: error term
# alpha0: intercept
#
# X: original design matrix (n_obs x p)
# X_nl: nonlinear terms design matrix (after Scheipl's decomposition & 
# reparameterization) (n_obs x \sum d_j)
# X_l: linear term design matrix (after Scheipl's decomposition & 
# reparameterization) (n_obs x p)
#
# beta_nl: vector that concatenates batch of coefficients for nonlinear 
# terms (\sum d_j x 1)
# beta_l: vector of coefficients for linear terms (p x 1)

# Indices 
# i: index of the observations in 1, ..., n_obs
# j: index of the covariates in 1, ..., p
# k, kk: instrumental (kk nested in k)

#' gendata: generative mechanism for simulating synthetic data
#' 
# @param n_obs number of observations (default to 200)
# @param p number of covariates (default to 10)
# @param sp proportion of nonzero linear terms, that is sparsity (default to .20)
# @param minb minimum effect for non null linear terms (default to 1.5)
# @param maxb maximum effect for non null linear terms (default to 3.0)
# @param error error term (default is 0.25)
# @return a list of 7 elements.
#   \itemize{
#   \item X - Design matrix (n_obs x p)
#   \item Xl - Matrix of linear terms (n_obs x p)
#   \item Xtilde - Matrix of linear terms (n_obs x  d*1)
#   \item d - vector of number of basis used to represent each covariate (p x 1)
#   \item intcp - scalat, value of the offset
#   \item betal - vector of coefficients for linear terms (p x 1)
#   \item betanl - vector that concatenates batch of coefficients for nonlinear terms (d*1 x 1)
# }
#' @export
gendata <- function(n_obs = 200, p = 10, minb = 1.5, maxb = 3.0, error = 0.01, 
                    scenario = 4, nnnc = 3, noi = 3, ha = 2) {
  # ha = 2, strong heredity assumption
  # ha = 1, weak heredity assumption
  # ha = 0, no assumption

  # noi non null interaction for each non null principal effect
  p_vector <- 1:p
  
  X <- vector()
  X_l <- vector()
  X_nl <- vector()
  d <- vector()
  
  for(j in 1:p) {
    xj <- rnorm(n_obs, 0, 1)
    X <- cbind(X, xj)
    xjl <- lin(xj)
    X_l <- cbind(X_l, xjl)
    xjtilde <- sm(x = xj, rankZ = .95) 
    X_nl <- cbind(X_nl, xjtilde)
    d[j] <- dim(xjtilde)[2]
  }
  cd <- c(0, cumsum(d))
  q <- sum(d)
  
  # Scenario 1
  # Linear main effect and linear interaction effect
  # Scenario 2
  # Linear main effect and non-linear interaction effect
  # Scenario 3
  # Non-linear main effect and linear interaction effect
  # Scenario 4
  # Non-linear main effect and non-linear interaction effect
  
  # alpha_0 main effect (linear)
  theta_l <- matrix(0, nrow = n_obs, ncol = p)
  # nnnc <- round(p*sp) # number of non null covariates (linear terms) (IMPORTANTE)
  innc <- sort(sample(1:p, nnnc, replace = FALSE)) # position of non null covariates (linear terms)
  for (j in innc) {
    theta_l[,j] <- runif(1, minb, maxb) * sign(runif(1, -1, 1))
  }
  
  # alpha_0 main effect (non linear)
  theta_nl <- matrix(0, nrow = n_obs, ncol = p)
  if ((scenario == 3) | (scenario == 4)) {
    innc_nl <- innc # indices for non null covariates (non linear terms)
    for (j in innc) {
      theta_nl[,j] <- runif(1, minb, maxb)*sign(runif(1, -1, 1))
      if (sign(theta_l[1,j]) != sign(theta_nl[1,j])) {
        theta_nl[,j] <- theta_nl[,j]*(-1)
      }
    }
  }
  # interaction omega parameter (beta star in the manuscript)
  omega_l <- matrix(0, p, p)
  # linnc <- length(innc)
  # interaction <- matrix(0, length(innc), 2)
  if (length(innc) != 1) {
    if (ha == 2) { # strong heredity
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          if ((theta_l[1,i] != 0) & (theta_l[1,j] != 0)) {
            omega_l[i,j] <- theta_l[1,i]*theta_l[1,j]
          }
        }
      }
    }
  }
  
  omega_nl <- matrix(0, nrow = p, ncol = q)
  if (ha == 1) { # weak heredity
    lastInt <- setdiff(p_vector, innc)
    for (i in innc) {
      innc_perm <- p_vector[p_vector != i]
      if (i == p) {
        nn_int <- sample(lastInt, 1)
      } else {
        nn_int <- sample(innc_perm, 1)
      }
      vectorInter <- sort(c(i,nn_int))
      omega_l[vectorInter[1],vectorInter[2]] <- rnorm(1, 2, 0.5)*sign(runif(1, -1, 1))
      if (scenario == 2 || scenario == 4) {
        omega_nl[vectorInter[1], (cd[vectorInter[2]]+1):(cd[vectorInter[2]+1])] <- rnorm(1, 2, 0.5)*sign(runif(1, -1, 1))
      }
    }
  }
  
  if (ha == 0) {
    cov1 <- sort(sample(1:p, noi, replace = FALSE))
    chCov2 <- setdiff(1:p, cov1)
    cov2 <- sort(sample(chCov2, noi, replace = TRUE))
    for (i in 1:noi) {
      first <- cov1[i]
      second <- cov2[i]
      vectorInter <- sort(c(first, second))
      first <- vectorInter[1]
      second <- vectorInter[2]
      omega_l[first, second] <- rnorm(1, 2, 0.5)*sign(runif(1, -1, 1))
      if (scenario == 2 || scenario == 4) {
        omega_nl[first, (cd[second]+1):(cd[second+1])] <- rnorm(1, 2, 0.5)*sign(runif(1, -1, 1))
      }
    }
  }
  
  # if (scenario == 4 || scenario == 2) {
  #   if (ha == 2) {
  #     for (i in 1:(p-1)) {
  #       for (j in (i+1):p) {
  #         if ((theta_nl[1,i] != 0) & (theta_nl[1,j] != 0)) {
  #           omega_nl[i, (cd[j]+1):(cd[j+1])] <- theta_nl[1,i]*theta_nl[1,j]
  #         }
  #       }
  #     }
  #   }
  # }
  
  if (ha == 2) {
    if (scenario == 2 || scenario == 4) {
      for (i in 1:(p-1)) {
        for (j in (i+1):p) {
          if ((theta_l[1,i] != 0) & (theta_l[1,j] != 0)) {
            omega_nl[i, (cd[j]+1):(cd[j+1])] <- rnorm(1, 2, 0.5)*sign(runif(1, -1, 1))
          }
        }
      }
    }
  }
  
  # alpha lin and non lin
  alpha_l <- matrix(0, nrow = n_obs, ncol = p)
  # compute alpha
  for (j in 1:p) {
    alpha_l[,j] <- theta_l[,j]
    if (j != p) {
      for (k in (j+1):p) {
        alpha_l[,j] <- alpha_l[,j] + matrix(X_l[,k])*omega_l[j,k]
      }
    }
  }
  
  alpha_nl <- matrix(0, nrow = n_obs, ncol = p)
  for (j in 1:p) {
    alpha_nl[,j] <- theta_nl[,j]
    if (j != p) {
      for (k in (j+1):p) {
        alpha_nl[,j] <- alpha_nl[,j] + X_nl[,(cd[k]+1):(cd[k+1])]%*%omega_nl[j, (cd[k]+1):(cd[k+1])]
      }
    }
  }
  
  # inizialize xi linear and xi non linear
  # xi_l <- sample(c(-1, 1), p, replace = TRUE)
  # xi_tilde <- c()
  # for (i in 1:p) {
  #   xi_tilde <- c(xi_tilde, rep(xi_l[i], d[i]))
  # }
  m_l <- sample(c(-1, 1), size = p, replace = TRUE, prob = c(0.5, 0.5))
  m_nl <- sample(c(-1, 1), size = q, replace = TRUE, prob = c(0.5, 0.5))
  xi_l <- rnorm(p, m_l, 1)
  xi_tilde <- rnorm(q, m_nl, 1)
  
  # compute beta linear and non linear
  beta_l <- matrix(0, nrow = n_obs, ncol = p)
  for (j in 1:p) {
    beta_l[,j] <- alpha_l[,j] * xi_l[j]
  }
  
  beta_nl <- matrix(0, nrow = n_obs, ncol = q)
  for (j in 1:p) {
    for(i in 1:n_obs) {
      beta_nl[i,(cd[j]+1):(cd[j+1])] <- alpha_nl[i,j]*xi_tilde[(cd[j]+1):(cd[j+1])]
    }
  }
  
  eta0 <- rep(runif(1, 1, 2)*sign(runif(1, -1, 1)), n_obs) #intercept
  epsilon <- rnorm(n_obs, 0, error) #error term
  
  Y <- rep(0, n_obs)
  for (i in 1:n_obs) {
    for (j in 1:p) {
      Y[i] <- Y[i] + X_l[i,j]*beta_l[i,j] + X_nl[i,(cd[j]+1):(cd[j+1])]%*%beta_nl[i,(cd[j]+1):(cd[j+1])]
    }
    Y[i] <- eta0[i] + Y[i] + epsilon[i]
  }

  return(list(
    Y = Y,
    X = X, 
    X_l = X_l, 
    X_nl = X_nl, 
    d = d, 
    intcpt = eta0[1],
    beta_l = beta_l, 
    beta_nl = beta_nl,
    theta_l = theta_l,
    theta_nl = theta_nl,
    alpha_l = alpha_l,
    alpha_nl = alpha_nl,
    omega_l = omega_l,
    omega_nl = omega_nl))
} # closes function genmech

#' my_indicies
my_indices <- function(est, truth, mpp = T){
  if(mpp == T){
    mppi <- apply(est, 2, mean)
    sel <- factor(as.numeric(mppi > 0.5), levels = c(0,1))
  } else {
    sel <- est
  }
  true <- factor(as.numeric(truth != 0), levels = c(0,1))
  contTable <- matrix(table(sel, true), nrow = 2, ncol = 2)
  matt <- mcc(confusionM = contTable)
  tpr <- contTable[1,1]/sum(contTable[,1]) 
  fpr <- contTable[1,2]/sum(contTable[,2]) 
  return(list(matt, tpr, fpr, sel))
}

#' my_indices_int
my_indices_int <- function(est, truth, linear = TRUE, d, omega_nl){
  p <- dim(est[,,1])[1]
  nout <- length(est)
  if (linear) {
    #est <- array(unlist(est), dim = c(p, p, nout))
    mppi <- apply(est, c(1,2), mean)
    sel <- mppi > 0.5
    sel[sel == TRUE] = 1
    true <-  truth != 0
    true[true == TRUE] = 1
  } else {
    q <- sum(d)
    cd <- c(0, cumsum(d))
    #est <- array(unlist(est), dim = c(p, p, nout))
    mppi <- apply(est, c(1,2), mean)
    sel <- mppi > 0.5
    sel[sel == TRUE] = 1
    TData <- omega_nl
    true <- matrix(0, nrow = p, ncol = p)
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        if (sum(TData[i, (cd[j]+1):(cd[j+1])]) != 0) {
          true[i,j] <- 1
        }
      } 
    }
  }
  
  tp <- 0
  tn <- 0
  fp <- 0
  fn <- 0
  
  for (j in 1:(p-1)) {
    for (jj in (j+1):p) {
      if ((sel[j,jj] == true[j,jj]) & (true[j,jj] == 1)) {
        tp <- tp + 1
      }
      if ((sel[j,jj] == true[j,jj]) & (true[j,jj] == 0)) {
        tn <- tn + 1
      }
      if ((sel[j,jj] != true[j,jj]) & (true[j,jj] == 0)) {
        fp <- fp + 1
      }
      if ((sel[j,jj] != true[j,jj]) & (true[j,jj] == 1)) {
        fn <- fn + 1
      }
    }
  }
  contTable <- matrix(c(tp, fp, fn, tn), nrow = 2, ncol = 2)
  matt <- mcc(confusionM = contTable)
  # True Positive Rate
  tpr <- contTable[1,1]/sum(contTable[,1])
  # False Positive Rate
  fpr <- contTable[1,2]/sum(contTable[,2])
  return(list(matt, tpr, fpr, sel))
}

######## PLOTS ###########
# Single MPPI Plot
#' @export
plot_mppi <- function(df, title) {
  # automatic step for the x-axis
  step <- ceiling(dim(df)[1]*0.20)
  selected_labels <- seq(from = 0, to = nrow(df), by = step)
  selected_labels[1] <- 1
  
  mainPlot <- ggplot(df, aes(x = cov, y = mppi)) +
    geom_segment(aes(xend = cov, yend = 0), color = 'black') +
    labs(x = 'Covariate Index', y = 'MPPI', title = title) +
    geom_hline(yintercept = 0.5, linetype = 'dashed', color = 'black') +
    theme_minimal() +
    theme(panel.grid.major = element_blank()) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_discrete(
      breaks = selected_labels,
      labels = selected_labels
    )
}

# Plot the MPPI
#' @export
mppi_plot <- function(resultMCMC) {
  # Function that returns the MPPI plot. The MCMC must be run with the option detail = TRUE.
  mainPlot <- NULL
  if (!is.null(resultMCMC$gamma_l)) {
    # stop("The MCMC should be run with detail = TRUE")
    p <- ncol(resultMCMC$gamma_star_l[[1]])
    nlp <- ncol(resultMCMC$gamma_star_nl[[1]])
    nout <- nrow(resultMCMC$linear_predictor)
    gammaStarLin <- array(unlist(resultMCMC$gamma_star_l), dim = c(p, p, nout))
    gammaStarNLin <- array(unlist(resultMCMC$gamma_star_nl), dim = c(nlp, nlp, nout))
    # gamma main linear
    gammaLin <- resultMCMC$gamma_l
    mppi_MainLinear <- apply(gammaLin, 2, mean)
    # gamma main non linear
    gammaNLin <- resultMCMC$gamma_nl
    mppi_MainNonLinear <- apply(gammaNLin, 2, mean)
    # gamma star linear (list of matrix)
    mppi_IntLinear <- apply(gammaStarLin, c(1,2), mean)
    # gamma star non linear
    mppi_IntNonLinear <- apply(gammaStarNLin, c(1,2), mean)
    # upper triangular matrix
    mppi_IntLinear <- mppi_IntLinear[upper.tri(mppi_IntLinear)]
    mppi_IntLinear <- mppi_IntLinear[mppi_IntLinear != 0]
    mppi_IntNonLinear <- mppi_IntNonLinear[upper.tri(mppi_IntNonLinear)]
    mainL <- data.frame(cov = as.factor(1:length(mppi_MainLinear)), mppi = mppi_MainLinear)
    mainNL <- data.frame(cov = as.factor(1:length(mppi_MainNonLinear)), mppi = mppi_MainNonLinear)
    interL <- data.frame(cov = as.factor(1:length(mppi_IntLinear)), mppi = mppi_IntLinear)
    interNL <- data.frame(cov = as.factor(1:length(mppi_IntNonLinear)), mppi = mppi_IntNonLinear)
    # Plot
    plot1 <- plot_mppi(mainL, 'Linear Main Effect')
    plot2 <- plot_mppi(mainNL, 'Non-Linear Main Effect')
    plot3 <- plot_mppi(interL, 'Linear Interaction Effect')
    plot4 <- plot_mppi(interNL, 'Non-Linear Interaction Effect')
    mainPlot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
  } else if (!is.null(resultMCMC$mppi_MainLinear)) {
    # Linear Interaction effects
    mppi_IntLinear <- resultMCMC$mppi_IntLinear
    mppi_IntLinear <- mppi_IntLinear[upper.tri(mppi_IntLinear)]
    # check null interaction 
    cNullInt <- mppi_IntLinear[mppi_IntLinear != 0]
    # NonLinear Interaction effects
    mppi_IntNonLinear <- resultMCMC$mppi_IntNonLinear
    mppi_IntNonLinear <- mppi_IntNonLinear[upper.tri(mppi_IntNonLinear)]
    mainL <- data.frame(cov = as.factor(1:length(resultMCMC$mppi_MainLinear)), mppi = resultMCMC$mppi_MainLinear)
    mainNL <- data.frame(cov = as.factor(1:length(resultMCMC$mppi_MainNonLinear)), mppi = resultMCMC$mppi_MainNonLinear)
    interL <- data.frame(cov = as.factor(1:length(cNullInt)), mppi = cNullInt)
    interNL <- data.frame(cov = as.factor(1:length(mppi_IntNonLinear)), mppi = mppi_IntNonLinear)
    # Plot
    plot1 <- plot_mppi(mainL, 'Linear Main Effect')
    plot2 <- plot_mppi(mainNL, 'Non-Linear Main Effect')
    plot3 <- plot_mppi(interL, 'Linear Interaction Effect')
    plot4 <- plot_mppi(interNL, 'Non-Linear Interaction Effect')
    mainPlot <- grid.arrange(plot1, plot2, plot3, plot4, nrow = 2, ncol = 2)
  } else {
    stop("I am unable to obtain any information about the MPPI.")
  }
}


# Main and Interaction Effects (Linear/Non-Linear)
#' @export
plotEffectResponse <- function(imod, idx, type = "linear", effect = "main",
                               xtitle = "", ytitle = "", ztitle = "", 
                               azimut = 30, polar = 20) {
  X_l <- imod$X_lin
  X_nl <- imod$X_nl
  ngf <- dim(X_nl)[1]
  d <- c(imod$d)
  pL <- dim(X_l)[2]
  pNL <- length(d)
  cd <- c(0, cumsum(d))
  q <- sum(d)
  nout <- length(imod$intercept)
  ones <- matrix(1, nrow = ngf, ncol = 1)
  resp <- matrix(0, nrow = ngf, ncol = length(imod$omega_l))
  if (length(idx) == 2) {
    idx1 <- idx[1]
    idx2 <- idx[2]
  }
  # Warnings
  if (length(idx) == 1 && effect == "main") {
    message("Plotting main effect for variable: ", idx)
    
  } else if (length(idx) == 2 && effect == "interaction") {
    message("Plotting interaction effect between: ", 
            paste0(idx1, " and ", idx2))
    
  } else {
    if (effect == "main") {
      stop("Main effect requires exactly 1 index. Received: ", length(idx))
    } else if (effect == "interaction") {
      stop("Interaction effect requires exactly 2 indices. Received: ", length(idx))
    } else {
      stop("Invalid effect type: ", effect, ". Must be 'main' or 'interaction'")
    }
  }
  # check MPPI
  if (effect == "main" & type == "linear") {
    cat("Linear Main Effects (MPPI >= 0.5):\n")
    mppi_linear <- colMeans(imod$gamma_l)
    sig_linear <- which(mppi_linear >= 0.5)
    if(length(sig_linear) > 0) {
      for(j in sig_linear) {
        cat("-", "Covariate", j, ":", round(mppi_linear[j], 2), "\n")
      }
    } else {
      stop("None found\n")
    }
  }
  
  if (effect == "main" & type == "nonlinear") {
    cat("\nNon-Linear Main Effects (MPPI >= 0.5):\n")
    mppi_nonlinear <- colMeans(imod$gamma_nl)
    sig_nonlinear <- which(mppi_nonlinear >= 0.5)
    if(length(sig_nonlinear) > 0) {
      for(j in sig_nonlinear) {
        cat("-", "Covariate", j, ":", round(mppi_nonlinear[j], 2), "\n")
      }
    } else {
      stop("None found\n")
    }
  }
  
  if (effect == "interaction" & type == "linear") {
    gammaStarLin <- array(unlist(imod$gamma_star_l), dim = c(pL, pL, nout))
    mppi_IntLinear <- apply(gammaStarLin, c(1,2), mean)
    var_names <- colnames(imod$X_lin) %||% paste0("V", 1:pL)
    cat("\nLinear Interaction Effects (MPPI >= 0.5):\n")
    sig_interactions <- which(mppi_IntLinear >= 0.5 & upper.tri(mppi_IntLinear), arr.ind = TRUE)
    if(nrow(sig_interactions) > 0) {
      for(i in 1:nrow(sig_interactions)) {
        row_idx <- sig_interactions[i, 1]
        col_idx <- sig_interactions[i, 2]
        cat("-", var_names[row_idx], "x", var_names[col_idx], 
            ":", round(mppi_IntLinear[row_idx, col_idx], 2), "\n")
      }
    } else {
      stop("No significant Linear interactions found\n")
    }

  }
  
  if (effect == "interaction" & type == "nonlinear") {
    gammaStarNLin <- array(unlist(imod$gamma_star_nl), dim = c(pNL, pNL, nout))
    mppi_IntNLinear <- apply(gammaStarNLin, c(1,2), mean)
    var_names <- colnames(imod$X_nl) %||% paste0("V", 1:pNL)
    cat("\nNon-Linear Interaction Effects (MPPI >= 0.5):\n")
    sig_interactions <- which(mppi_IntNLinear >= 0.5 & upper.tri(mppi_IntNLinear), arr.ind = TRUE)
    if(nrow(sig_interactions) > 0) {
      for(i in 1:nrow(sig_interactions)) {
        row_idx <- sig_interactions[i, 1]
        col_idx <- sig_interactions[i, 2]
        cat("-", var_names[row_idx], "x", var_names[col_idx], 
            ":", round(mppi_IntNLinear[row_idx, col_idx], 2), "\n")
      }
    } else {
      stop("No significant Non-Linear interactions found\n")
    }
  }
  
  for (t in 1:length(imod$omega_l)) {
    # omega linear and nonlinear
    omega_l <- imod$omega_l[[t]]
    omega_nl <- imod$omega_nl[[t]]
    # main linear and non linear effects
    theta_l <- imod$theta_l[t,]
    theta_nl <- imod$theta_nl[t,]
    # inizialize xi linear and xi non linear
    xi_l <- imod$xi_l[t,]
    xi_nl <- imod$xi_nl[t,]
    resp[,t] <- imod$intercept[t]*ones
    # Main effect Linear
    for (j in 1:pL) {
      if (type == "linear" & effect == "main") {
        if (j == idx) {
          resp[,t] <- resp[,t] + matrix(X_l[,j])*theta_l[j]*xi_l[j]
        } else {
          resp[,t] <- resp[,t] + median(matrix(X_l[,j])*theta_l[j]*xi_l[j])
        }
      } else {
        resp[,t] <- resp[,t] + median(matrix(X_l[,j])*theta_l[j]*xi_l[j])
      }
    }
    # Main effect NonLinear
    for (j in 1:pNL) {
      if (type == "nonlinear" & effect == "main") {
        if (j == idx) {
          resp[,t] <- resp[,t] + rowSums(X_nl[,(cd[j]+1):(cd[j+1])]*(as.matrix(theta_nl[j]*ones)%*%as.matrix(t(xi_nl[(cd[j]+1):(cd[j+1])]))))
        } else {
          resp[,t] <- resp[,t] + median(rowSums(X_nl[,(cd[j]+1):(cd[j+1])]*(as.matrix(theta_nl[j]*ones)%*%as.matrix(t(xi_nl[(cd[j]+1):(cd[j+1])])))))
        }
      } else {
        resp[,t] <- resp[,t] + median(rowSums(X_nl[,(cd[j]+1):(cd[j+1])]*(as.matrix(theta_nl[j]*ones)%*%as.matrix(t(xi_nl[(cd[j]+1):(cd[j+1])])))))
      }
    }
    # Interaction effect Linear
    for (j in 1:(pL-1)) {
      for (k in (j+1):pL) {
        # linear
        if (type == "linear" & effect == "interaction") {
          if (j == idx1 & k == idx2) {
            resp[,t] <- resp[,t] + matrix(X_l[,j])*matrix(X_l[,k])*omega_l[j,k]*xi_l[j]
          } else {
            resp[,t] <- resp[,t] + median(matrix(X_l[,j])*matrix(X_l[,k])*omega_l[j,k]*xi_l[j])
          }
        } else {
          resp[,t] <- resp[,t] + median(matrix(X_l[,j])*matrix(X_l[,k])*omega_l[j,k]*xi_l[j])
        }
      }
    }
    # Interaction effect NonLinear
    for (j in 1:(pNL-1)) {
      for (k in (j+1):pNL) {
        # nonlinear
        if (type == "nonlinear" & effect == "interaction") {
          if (j == idx1 & k == idx2) {
            resp[,t] <- resp[,t] + rowSums(X_nl[,(cd[j]+1):(cd[j+1])]*(X_nl[,(cd[k]+1):(cd[k+1])]%*%as.matrix(omega_nl[j, (cd[k]+1):(cd[k+1])]))%*%as.matrix(t(xi_nl[(cd[j]+1):(cd[j+1])])))
          } else {
            resp[,t] <- resp[,t] + median(rowSums(X_nl[,(cd[j]+1):(cd[j+1])]*(X_nl[,(cd[k]+1):(cd[k+1])]%*%as.matrix(omega_nl[j, (cd[k]+1):(cd[k+1])]))%*%as.matrix(t(xi_nl[(cd[j]+1):(cd[j+1])]))))
          }
        } else {
          resp[,t] <- resp[,t] + median(rowSums(X_nl[,(cd[j]+1):(cd[j+1])]*(X_nl[,(cd[k]+1):(cd[k+1])]%*%as.matrix(omega_nl[j, (cd[k]+1):(cd[k+1])]))%*%as.matrix(t(xi_nl[(cd[j]+1):(cd[j+1])]))))
        }
      }
    }
  }
  # PLOT
  if (effect == "main") {
    xj <- X_l[, idx]
    y_median <- apply(resp, 1, median)
    y_low <- apply(resp, 1, function(x) quantile(x, 0.025))
    y_up <- apply(resp, 1, function(x) quantile(x, 0.975))
    df_smooth <- data.frame(
      xj = xj,
      yx = predict(smooth.spline(xj, y_median, spar = 0.7), x = xj)$y,
      yx_low = predict(smooth.spline(xj, y_low, spar = 0.7), x = xj)$y,
      yx_up = predict(smooth.spline(xj, y_up, spar = 0.7), x = xj)$y
    )
    plotRend <- ggplot(df_smooth, aes(x = xj)) +
      geom_ribbon(aes(ymin = y_low, ymax = y_up), 
                  fill = "black", 
                  alpha = 0.2,   
                  color = NA) +  
      geom_line(aes(y = yx),
                color = "black", 
                linewidth = 0.8,     
                alpha = 0.9) +  
      geom_rug(
        data = data.frame(x_observed = xj),
        aes(x = x_observed),
        sides = "b",
        color = "black",  
        alpha = 0.3,
        inherit.aes = FALSE
      ) +
      ylab(ytitle) + 
      xlab(xtitle) + 
      theme_minimal() +
      theme(
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_blank()
      )
    print(plotRend) 
  } else {
    y_median <- apply(resp, 1, median)
    y_low <- apply(resp, 1, function(x) quantile(x, 0.025))
    y_up <- apply(resp, 1, function(x) quantile(x, 0.975))
    xj <- X_l[, idx1]
    xk <- X_l[, idx2]
    nx <- 200
    ny <- 200
    xo <- seq(min(xj), max(xj), length.out = nx)
    yo <- seq(min(xk), max(xk), length.out = ny)
    interp_safe <- function(x, y, z) {
      tryCatch({
        interp(x = x, y = y, z = z, 
               xo = xo, yo = yo,
               duplicate = "mean",
               linear = FALSE,
               extrap = FALSE)$z
      }, error = function(e) {
        matrix(NA, nrow = nx, ncol = ny)
      })
    }
    z_median <- interp_safe(xj, xk, y_median)
    z_low <- interp_safe(xj, xk, y_low)
    z_up <- interp_safe(xj, xk, y_up)
    stopifnot(
      all(dim(z_median) == c(nx, ny)),
      all(dim(z_low) == c(nx, ny)),
      all(dim(z_up) == c(nx, ny))
    )
    grid <- mesh(xo, yo)
    # Okabe-Ito
    okabe_ito_continuous <- colorRampPalette(c("#F0E442", "#E69F00", 
                                               "#009E73", "#D55E00", 
                                               "#CC79A7", "#0072B2", 
                                               "#999999", "#000000"))
    col_palette <- okabe_ito_continuous(100)
    grid_color <- "grey90"
    par(mar = c(2, 2, 2, 2),
        cex.axis = 0.8,     
        cex.lab = 0.9)
    surf3D(
      x = grid$x, 
      y = grid$y,
      z = z_median,
      col = col_palette,
      colkey = FALSE, 
      lighting = list(ambient = 0.3, diffuse = 0.7, specular = 0.4),
      shade = 0.1,
      theta = azimut,
      phi = polar,
      bty = "b2",
      ticktype = "detailed",
      xlab = xtitle,
      ylab = ytitle,
      zlab = ztitle,
      zlim = range(z_median, na.rm = TRUE),
      resfac = 1,
      border = NA,
      facets = TRUE,
      col.grid = grid_color,
      axes = TRUE,
      box = TRUE,
      rasterImage = TRUE,
      cex.axis = 0.6,
      cex.lab = 1
    )
    contour3D(
      x = xo,
      y = yo,
      z = z_median,
      colvar = z_median,
      add = TRUE,                   
      col = "black",                
      lwd = 0.8,                    
      nlevels = 10,                 
      drawlabels = TRUE,            
      labcex = 0.7,                 
      alpha = 0.7,                  
      method = "flattest"           
    )
  }
  if (effect == "main") {
    return(invisible(list(
      plot = plotRend,
      data = df_smooth,
      variable = idx,
      effect_type = effect
    )))
  } else {
    return(invisible(list(
      surface_data = list(
        x = grid$x,
        y = grid$y,
        z_median = z_median,
        z_low = z_low,
        z_up = z_up
      ),
      variables = c(idx1, idx2),
      effect_type = effect
    )))
  }
}


# Covariate effect on Beta regression coefficient (Linear/Non-Linear)
#' @export
plotEffectBeta <- function(imod, idx, modifier, type = "linear",
                           xtitle = "", ytitle = "") {
  X_l <- imod$X_lin
  X_nl <- imod$X_nl
  ngf <- dim(X_nl)[1]
  d <- c(imod$d)
  pL <- dim(X_l)[2]
  pNL <- length(d) 
  cd <- c(0, cumsum(d))
  q <- sum(d)
  ones <- matrix(1, nrow = ngf, ncol = 1)
  Beta_LinearIdx <- matrix(0, nrow = ngf, ncol = length(imod$omega_l))
  Beta_NonLinearIdx <- matrix(0, nrow = ngf, ncol = length(imod$omega_l))
  for (t in 1:length(imod$omega_l)) {
    # omega linear and nonlinear
    omega_l <- imod$omega_l[[t]]
    omega_nl <- imod$omega_nl[[t]]
    # main linear and non linear effects
    theta_l <- imod$theta_l[t,]
    theta_nl <- imod$theta_nl[t,]
    # inizialize xi linear and xi non linear
    xi_l <- imod$xi_l[t,]
    xi_nl <- imod$xi_nl[t,]
    if (type == "linear") {
      Beta_LinearIdx[,t] <- theta_l[idx]*ones*xi_l[idx]
      if (idx != pL) {
        for (k in (idx+1):pL) {
          if (k == modifier) {
            Beta_LinearIdx[,t] <- Beta_LinearIdx[,t] + X_l[,k]*omega_l[idx,k]*xi_l[idx]
          } else {
            Beta_LinearIdx[,t] <- Beta_LinearIdx[,t] + median(X_l[,k]*omega_l[idx,k]*xi_l[idx])
          }
        }
      }
    } else {
      Beta_NonLinearIdx[,t] <- median(rowSums(theta_nl[idx]*ones%*%as.matrix(t(xi_nl[(cd[idx]+1):(cd[idx+1])]))))
      if (idx != pNL) {
        for (k in (idx+1):pNL) {
          if (k == modifier) {
            Beta_NonLinearIdx[,t] <- Beta_NonLinearIdx[,t] + rowSums(X_nl[,(cd[k]+1):(cd[k+1])]%*%as.matrix(omega_nl[idx, (cd[k]+1):(cd[k+1])])%*%as.matrix(t(xi_nl[(cd[idx]+1):(cd[idx+1])])))
          } else {
            Beta_NonLinearIdx[,t] <- Beta_NonLinearIdx[,t] + median(rowSums(X_nl[,(cd[k]+1):(cd[k+1])]%*%as.matrix(omega_nl[idx, (cd[k]+1):(cd[k+1])])%*%as.matrix(t(xi_nl[(cd[idx]+1):(cd[idx+1])]))))
          }
        } 
      }
    }
  }
  xj <- X_l[, modifier]
  if (type == "linear") {
    betaLin_median <- apply(Beta_LinearIdx, 1, median)
    betaLin_low <- apply(Beta_LinearIdx, 1, function(x) quantile(x, 0.025))
    betaLin_up <- apply(Beta_LinearIdx, 1, function(x) quantile(x, 0.975))
    df_smooth <- data.frame(
      xj = xj,
      beta = predict(smooth.spline(xj, betaLin_median, spar = 0.7), x = xj)$y,
      beta_low = predict(smooth.spline(xj, betaLin_low, spar = 0.7), x = xj)$y,
      beta_up = predict(smooth.spline(xj, betaLin_up, spar = 0.7), x = xj)$y
    )
  } else {
    betaNLin_median <- apply(Beta_NonLinearIdx, 1, median)
    betaNLin_low <- apply(Beta_NonLinearIdx, 1, function(x) quantile(x, 0.025))
    betaNLin_up <- apply(Beta_NonLinearIdx, 1, function(x) quantile(x, 0.975))
    df_smooth <- data.frame(
      xj = xj,
      beta = predict(smooth.spline(xj, betaNLin_median, spar = 0.7), x = xj)$y,
      beta_low = predict(smooth.spline(xj, betaNLin_low, spar = 0.7), x = xj)$y,
      beta_up = predict(smooth.spline(xj, betaNLin_up, spar = 0.7), x = xj)$y
    )
  }
  plotRend <- ggplot(df_smooth, aes(x = xj)) +
    geom_ribbon(aes(ymin = beta_low, ymax = beta_up), 
                fill = "black", 
                alpha = 0.2,   
                color = NA) +  
    geom_line(aes(y = beta),
              color = "black", 
              linewidth = 0.8,     
              alpha = 0.9) +  
    geom_rug(
      data = data.frame(x_observed = xj),
      aes(x = x_observed),
      sides = "b",
      color = "black",  
      alpha = 0.3,
      inherit.aes = FALSE
    ) +
    ylab(ytitle) + 
    xlab(xtitle) + 
    theme_minimal() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    )
  print(plotRend) 
  return(plotRend)
}

