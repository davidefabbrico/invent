#' invMCMC
#' 
#' @export

##### ------------------------------------------------------------------ ######
invMCMC <- function(y, x, hyperpar = c(5, 25, 5, 5, 0.00025, 0.4, 1.6, 0.2, 1.8), 
                   mht = c(0.8, 1.4, 0.2, 1.4, 0.7, 0.15, 0.8, 0.4), 
                   rank = 0.95, iter = 10000, burnin = iter/2, thin = 5, ha = 2) {

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
  
  result = bodyMCMC(as.vector(y), as.integer(p), as.integer(nobs), as.vector(cd), 
                    as.vector(d), as.matrix(X_l), as.matrix(X_nl), 
                    as.vector(hyperpar), as.vector(mht), as.integer(iter), 
                    as.integer(burnin), as.integer(thin), ha)
  
  res <- result
  
  return(res)
}
