rm(list=ls())
set.seed(221)
devtools::load_all()

data_list <- gendata(n_obs = 50, p = 10, scenario = 1, nnnc = 3)
out <- nlivcC(data_list$Y, data_list$X, iter = 100, burnin = 50, thin = 1)