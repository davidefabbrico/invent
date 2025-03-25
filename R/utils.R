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
  step <- floor(dim(df)[1]*0.20)
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


# Non-Linear Interaction Plot


