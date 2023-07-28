#' invMCMC
#' 
#' @export

##### ------------------------------------------------------------------ ######
invMCMC <- function(y, x, x_val, hyperpar = c(5, 25, 5, 5, 0.00025, 0.4, 1.6, 0.2, 1.8), 
                   mht = c(1.4, 0.8, 1, 0.3, 0.7, 0.4, 4, 2.5), 
                   rank = 0.95, iter = 10000, burnin = iter/2, thin = 5, ha = 2, n_val=100) {

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
  pred <- TRUE
  if (n_val != 0) {
    x_val <- cbind(X_l[(nobs-n_val+1):nobs,], X_nl[(nobs-n_val+1):nobs,])
    X_l <- X_l[1:(nobs-n_val),]
    X_nl <- X_nl[1:(nobs-n_val),]
  } else {
    x_val <- matrix(0, nrow = 1, ncol = 1)
    pred <- FALSE
  }
  
  result = bodyMCMC(as.vector(y), as.integer(p), as.integer(nobs), as.vector(cd), 
                    as.vector(d), as.matrix(X_l), as.matrix(X_nl), 
                    as.vector(hyperpar), as.vector(mht), as.integer(iter), 
                    as.integer(burnin), as.integer(thin), ha, as.matrix(x_val), as.integer(pred))
  
  res <- result
  
  return(res)
}
