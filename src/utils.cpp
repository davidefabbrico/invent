#include <RcppArmadillo.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <stdlib.h> 
#include <chrono>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec myRange(int start, int end) {
  int dimension = end-start+1;
  arma::vec myVec(dimension);
  int j = 0;
  for(int i=start; i<=end; i++) {
    myVec(j) = i;
    j++;
  }
  return myVec;
}

// [[Rcpp::export]]
NumericVector callrgamma(int n, double shape, double scale) { 
  return(rgamma(n, shape, scale)); 
}

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}

// INTERCEPT
// [[Rcpp::export]]
double updateInterceptC(arma::vec y, int nobs, arma::vec lp_noInt, double sigma) {
  double intercept = R::rnorm(arma::accu(y-lp_noInt)/(nobs+1), sigma/(nobs+1));
  return intercept;
}

// M scalar
// [[Rcpp::export]]
double update_mCsca(double xi) {
  double m = 1/(1+exp(-2*xi));
  return m;
}

// M vector
// [[Rcpp::export]]
arma::vec update_mCvec(arma::rowvec xi) {
  int n = arma::size(xi)(1);
  arma::vec m(n);
  for (int i = 0; i<n; i++) {
    m(i) = 1/(1+exp(-2*xi(i)));
  }
  return m;
}

// TAU
// [[Rcpp::export]]
double update_tauC(double at, double bt, double alpha, double gamma) {
  double tau = 1/(callrgamma(1, at+0.5, 1/(bt + pow(alpha, 2)/(2*gamma)))(0));
  return tau;
}

//DELTA F Non sctructured
// [[Rcpp::export]]
int deltaFNSC(arma::vec x, double val) {
  int count = 0;
  int n = x.n_elem;
  for (int i = 0; i < n; i++) {
    if (x[i] == val) {
      count = count + 1;
    }
  }
  return count;
}

//DELTA F
// [[Rcpp::export]]
int deltaFSC(arma::mat x, double val) {
  int count = 0;
  int n = arma::size(x)(0);
  for (int j = 0; j<(n-1); j++) {
    for (int k = (j+1); k<n;k++) {
      if (x(j,k) == val) {
        count = count + 1;
      }
    }
  }
  return count;
}

// PI Non structured
// [[Rcpp::export]]
double update_piNSC(double ap, double bp, arma::vec gamma, double v0) {
  double pi = R::rbeta(ap + deltaFNSC(gamma, 1), bp + deltaFNSC(gamma, v0));
  return pi;
}

// PI Structured
// [[Rcpp::export]]
double update_piSC(double ap, double bp, arma::mat gamma, double v0) {
  double pi = R::rbeta(ap + deltaFSC(gamma, 1), bp + deltaFSC(gamma, v0));
  return pi;
}

// SIGMA
// [[Rcpp::export]]
double update_sigmaC(arma::colvec y, arma::colvec eta_pl, double as, double bs, int nobs) {
  double sigma = 1/(callrgamma(1, as + nobs/2, 1/(bs + arma::accu(pow((y - eta_pl), 2))/2))(0));
  return sigma;
}

// GAMMA Scalar
// [[Rcpp::export]]
double update_gammaScaC(double pi, double v0, double alpha, double tau) {
  double gamma_out = 0;
  double prob = 0;
  double d1 = log(pi/(1-pi)) + log(v0) * 1/2;
  double d2 = (1 - v0) / (2 * v0);
  double out = exp(d1 + d2 * pow(alpha, 2)/tau);
  if (out > 10000) {
    prob = 1;
    gamma_out = 1;
  } else {
    if (out < 0.00001) {
      prob = 0;
      gamma_out = v0;
    } else {
      prob = out / (1 + out);
      double samp = R::runif(0,1);
      if (samp < prob) {
        gamma_out = 1;
      } else {
        gamma_out = v0;
      }
    }
  }
  return gamma_out;
}

// GAMMA Vector
// [[Rcpp::export]]
arma::vec update_gammaVecC(double pi, double v0, arma::vec alpha, arma::vec tau) {
  int p = alpha.n_elem;
  arma::vec gamma_out(p);
  double prob = 0;
  double d1 = log(pi/(1-pi)) + log(v0) * 1/2;
  double d2 = (1 - v0) / (2 * v0);
  for (int j = 0; j<p; j++) {
    double out = exp(d1 + d2 * pow(alpha(j), 2)/tau(j));
    if (out > 10000) {
      prob = 1;
      gamma_out(j) = 1;
    } else {
      if (out < 0.00001) {
        prob = 0;
        gamma_out(j) = v0;
      } else {
        prob = out / (1 + out);
        double samp = R::runif(0,1);
        if (samp < prob) {
          gamma_out(j) = 1;
        } else {
          gamma_out(j) = v0;
        }
      }
    }
  }
  return gamma_out;
}

// [[Rcpp::export]]
arma::vec dnormLogVec(arma::vec x, arma::vec means, double sds) {
  int n = x.n_elem;
  arma::vec res(n);
  for(int i = 0; i < n; i++) {
    res(i) = R::dnorm(x(i), means(i), sds, TRUE);
  }
  return res;
}

// ALPHA
// [[Rcpp::export]]
List update_alphaC(arma::vec y, double sigma, double tau, double gamma, arma::vec eta_star, arma::vec eta_pl, double alpha_star, double alpha) {
  int acc_alpha = 0;
  double a = arma::accu(dnormLogVec(y, eta_star, sqrt(sigma)));
  double b = R::dnorm(alpha_star, 0, sqrt(tau*gamma), TRUE);
  double num_alpha = a+b;
  double den_alpha = arma::accu(dnormLogVec(y, eta_pl, sqrt(sigma))) + R::dnorm(alpha, 0, sqrt(tau*gamma), TRUE);
  double ratio_alpha = num_alpha - den_alpha;
  double lsamp = log(R::runif(0, 1));
  double alpha_new = 0;
  if (ratio_alpha > lsamp) {
    acc_alpha = acc_alpha + 1;
    alpha_new = alpha_star;
  } else {
    alpha_new = alpha;
  }
  return List::create(Named("alpha") = alpha_new,
                      Named("acceptance") = acc_alpha);
}

// XI Linear
// [[Rcpp::export]]
List update_xiLC(arma::vec y, arma::colvec eta_star, arma::colvec eta_pl, double sigma, double m, double xi_star, double xi) {
  int acc_xi = 0;
  double num_xi = arma::accu(dnormLogVec(y, eta_star, sqrt(sigma))) + R::dnorm(xi_star, m, 1, TRUE);
  double den_xi = arma::accu(dnormLogVec(y, eta_pl, sqrt(sigma))) + R::dnorm(xi, m, 1, TRUE);
  double ratio_xi = num_xi - den_xi;
  double lsamp = log(R::runif(0, 1));
  double xi_new = 0;
  if (ratio_xi > lsamp) {
    acc_xi = acc_xi + 1;
    xi_new = xi_star;
  } else {
    xi_new = xi;
  }
  
  return List::create(Named("xi") = xi_new,
                      Named("acceptance") = acc_xi);
}

// XI NonLinear
// [[Rcpp::export]]
List update_xiNLC(arma::vec y, arma::vec eta_star, arma::vec eta_pl, double sigma, arma::vec m, arma::vec xi_star, arma::vec xi) {
  int dj = xi_star.n_elem;
  arma::vec acc_xi(dj);
  double num_xi = arma::accu(dnormLogVec(y, eta_star, sqrt(sigma))) + arma::accu(dnormLogVec(xi_star, m, 1));
  double den_xi = arma::accu(dnormLogVec(y, eta_pl, sqrt(sigma))) + arma::accu(dnormLogVec(xi, m, 1));
  double ratio_xi = num_xi - den_xi;
  double lsamp = log(R::runif(0, 1));
  arma::vec xi_new;
  if (ratio_xi > lsamp) {
    acc_xi = acc_xi + 1;
    xi_new = xi_star;
  } else {
    xi_new = xi;
  }
  
  return List::create(Named("xi") = xi_new,
                      Named("acceptance") = acc_xi);
}

// sign function
// [[Rcpp::export]]
int mysign(double x) {
  if (x < 0) {
    return -1;
  } else {
    return 1;
  }
}

// Compute Linear Predictor
// [[Rcpp::export]]
arma::vec compLinPred(int nobs, int p, arma::vec cd, double eta0, arma::mat X_l, arma::mat beta_l, arma::mat X_nl, arma::mat beta_nl) {
  arma::vec Eta0(nobs);
  arma::vec eta_pl = Eta0.fill(eta0);
  for (int i = 0; i<nobs; i++) {
    for (int j = 0; j<p; j++) {
      eta_pl(i) = eta_pl(i) + X_l(i,j)*beta_l(i,j) + as_scalar(X_nl(i, span(cd[j], cd[j+1]-1))*beta_nl(i, span(cd[j], cd[j+1]-1)).t());
    }
  }
  return eta_pl;
}

// Body MCMC
// [[Rcpp::export]]
List bodyMCMC(arma::vec y, int p, int nobs, arma::vec cd, arma::vec d, arma::mat X_l, arma::mat X_nl, arma::vec hyperpar, arma::vec mht, int iter, int burnin, int thin, int ha) {
  // Time 
  // auto start = std::chrono::high_resolution_clock::now();
  ////////////////////////////////////////////////////
  ////////////////// Initial value //////////////////
  ///////////////////////////////////////////////////
  // no assumption = 0
  // wh = 1
  // sh = 2
  int q = arma::accu(d);
  arma::vec vecOnes = ones(nobs, 1);
  // pi star linear 
  double pi_star_l = R::rbeta(1,1);
  // pi star non linear
  double pi_star_nl = R::rbeta(1,1);
  // pi 0 linear
  double pi_0_l = R::rbeta(1,1);
  // pi 0 non linear
  double pi_0_nl = R::rbeta(1,1);
  // gamma star linear
  // upper trinagular matrix
  arma::mat gamma_star_l(p,p);
  for (int i = 0; i<(p-1); i++) {
    for (int j = (i+1); j<p; j++) {
      int randNumber = round(R::runif(0,1));
      gamma_star_l(i,j) = randNumber;
      if (randNumber == 0) {
        gamma_star_l(i,j) = hyperpar(4);
      }
    }
  }
  // gamma star non linear
  // upper triangular matrix 
  arma::mat gamma_star_nl = gamma_star_l;
  // gamma 0 linear (vector of dimension p)
  arma::vec gamma_0_l(p);
  for (int j = 0; j<p; j++) {
    int randNumber = round(R::runif(0,1));
    gamma_0_l(j) = hyperpar(4);
    if (randNumber == 0) {
      gamma_0_l(j) = randNumber;
    }
  }
  // gamma non linear (vector of dimension p)
  arma::vec gamma_0_nl = gamma_0_l;
  // tau star linear
  arma::mat tau_star_l(p,p);
  for (int i = 0; i<p-1; i++) {
    for (int j = i+1; j<p; j++) {
      tau_star_l(i,j) = 1/callrgamma(1,1,1)(0);
    }
  }
  // tau star non linear
  arma::mat tau_star_nl = tau_star_l;
  // tau linear (variance vector)
  arma::vec tau_0_l(p);
  for (int j = 0; j<p; j++) {
    tau_0_l(j) = 1/callrgamma(1, 1, 1)(0);
  }
  // tau non linear (variance vector)
  arma::vec tau_0_nl = tau_0_l;
  // alpha star linear square matrix of dimension p (for interaction)
  arma::mat alpha_star_l(p,p);
  for (int i = 0; i<p-1; i++) {
    for (int j = i+1; j<p; j++) {
      alpha_star_l(i,j) = R::rnorm(0,1);
    }
  }
  arma::mat alpha_star_nl = alpha_star_l;
  // xi star linear (square upper trinangular matrix of dimension p)
  arma::mat xi_star_l(p,p);
  for (int i = 0; i<(p-1); i++) {
    for (int j = i+1; j<p; j++) {
      xi_star_l(i,j) = 1;
    }
  }
  // xi star non linear, matrix of dimension pxq where q is the sum of the basis
  arma::mat xi_star_nl(p,q, fill::ones);
  // omega linear
  arma::mat omega_l(p,p);
  // omega non linear
  arma::mat omega_nl(p,q);
  // compute the omega linear and non linear
  for (int j = 0; j<(p-1); j++) {
    for (int k = (j+1); k<p; k++) {
      omega_l(j,k) = xi_star_l(j,k) * alpha_star_l(j,k);
      // omega non linear
      omega_nl(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1))*alpha_star_nl(j,k);
    }
  }
  // alpha 0 linear vector (theta in thesis)
  arma::vec alpha_0_l = as<arma::vec>(wrap(Rcpp::rnorm(p, 0, 1)));
  // alpha 0 non linear vector (theta in thesis)
  arma::vec alpha_0_nl = alpha_0_l;
  // alpha linear
  arma::mat alpha_l(nobs, p);
  // alpha non linear
  arma::mat alpha_nl(nobs, p);
  // compute the alpha linear and non linear element
  for (int j = 0; j<p; j++) {
    alpha_l.col(j) = alpha_0_l(j)*vecOnes;
    alpha_nl.col(j) = alpha_0_nl(j)*vecOnes;
    if (j != (p-1)) {
      for (int k = (j+1); k<p; k++) {
        alpha_l.col(j) = alpha_l.col(j) + X_l.col(j)*omega_l(j,k);
        // matrix multiplication and slice
        arma::vec col_ind = myRange(cd[k], cd[k+1]-1);
        arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
        alpha_nl.col(j) = alpha_nl.col(j) + X_nl.cols(ucol_ind)*omega_nl(j, span(cd[k], cd[k+1]-1)).t();
      }
    }
  }
  // m star linear matrix of dimension pxp
  arma::mat m_star_l(p,p);
  for (int j = 0; j<p; j++) {
    for (int k = (j+1); k<p; k++) {
      m_star_l(j,k) = mysign(R::runif(-1,1));
    }
  }
  // m star non linear 
  arma::mat m_star_nl(p,q);
  for (int j = 0; j<p; j++) {
    for (int k = 0; k<q; k++) {
      m_star_nl(j,k) = mysign(R::runif(-1,1));
    }
  }
  // matrix with sign of dimension pxq
  // m linear vector
  arma::vec m_l(p);
  for (int j = 0; j<p; j++) {
    m_l(j) = mysign(R::runif(-1, 1));
  }
  // m non linear vector
  arma::vec m_nl(q);
  for (int j = 0; j<q; j++) {
    m_nl(j) = mysign(R::runif(-1, 1));
  }
  // xi linear
  arma::vec xi_l = as<arma::vec>(wrap(Rcpp::rnorm(p, 0, 1)));
  // xi non linear
  arma::vec xi_nl = as<arma::vec>(wrap(Rcpp::rnorm(q, 0, 1)));
  // beta linear
  arma::mat beta_l(nobs, p);
  for (int j = 0; j<p; j++) {
    beta_l.col(j) = xi_l(j)*alpha_l.col(j);
  }
  // beta non linear
  arma::mat beta_nl(nobs, q);
  // iter also on the number of observations
  for (int j = 0; j<p; j++) {
    for (int i = 0; i<nobs; i++) {
      beta_nl(i, span(cd[j], cd[j+1]-1)) = xi_nl(span(cd[j], cd[j+1]-1)).t()*alpha_nl(i,j);
    }
  }
  // intercept
  double eta0 = R::rnorm(0, 1);
  // variance of the normal model
  double sigma = 1/callrgamma(1, 1, 1)(0);
  // Linear Predictor
  arma::vec eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
  // index for storage
  int nout = (iter - burnin)/thin;
  int idx = 0;
  ////////////////////////////////////////////////////
  //////////////// End Initial value ////////////////
  ///////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////
  /////////////// Start Store Results ////////////////
  ///////////////////////////////////////////////////
  // alpha star acc
  arma::mat alpha_star_l_acc(p,p);
  arma::mat alpha_star_nl_acc(p,p);
  // xi star acc
  arma::mat xi_star_l_acc(p,p);
  arma::mat xi_star_nl_acc(p,q);
  // alpha 0 (theta in thesis) acc
  arma::vec alpha_0_l_acc(p);
  arma::vec alpha_0_nl_acc(p);
  // xi acc
  arma::vec xi_l_acc(p);
  arma::vec xi_nl_acc(q);
  // pi star linear and non linear
  arma::vec PI_S_l(nout);
  arma::vec PI_S_nl(nout);
  // pi 0 linear and non linear
  arma::vec PI_0_l(nout);
  arma::vec PI_0_nl(nout);
  // linear predictor
  arma::mat ETA_PL(nout, nobs);
  // intercept
  arma::vec ETA0(nout);
  // sigma variance model
  arma::vec SIGMA(nout);
  // log-likelihood
  arma::vec LOGLIKELIHOOD(nout);
  // Gamma star linear and non linear
  List GAMMA_S_l(nout);
  List GAMMA_S_nl(nout);
  // Tau star linear and non linear
  List TAU_S_l(nout);
  List TAU_S_nl(nout);
  // Alpha star linear and non linear
  List ALPHA_S_l(nout);
  List ALPHA_S_nl(nout);
  // m star linear and non linear
  List M_S_l(nout);
  List M_S_nl(nout);
  // Xi star linear and non linear
  List XI_S_l(nout);
  List XI_S_nl(nout);
  // Omega linear and non linear
  List OMEGA_l(nout);
  List OMEGA_nl(nout);
  // m linear and non linear
  arma::mat M_l(nout, p);
  arma::mat M_nl(nout, q);
  // xi linear and non linear
  arma::mat XI_l(nout, p);
  arma::mat XI_nl(nout, q);
  // Beta linear and non linear
  List BETA_l(nout);
  List BETA_nl(nout);
  // Gamma star linear and non linear
  arma::mat GAMMA_0_l(nout, p);
  arma::mat GAMMA_0_nl(nout, p);
  // Tau 0 linear and non linear
  arma::mat TAU_0_l(nout, p);
  arma::mat TAU_0_nl(nout, p);
  // Alpha 0 linear and non linear
  arma::mat ALPHA_0_l(nout, p);
  arma::mat ALPHA_0_nl(nout, p);
  // Alpha linear and non linear
  List ALPHA_l(nout);
  List ALPHA_nl(nout);
  // Profiling Matrix
  arma::mat Prof_Mat(nout, 2);
  // Init Parameters
  double alpha_star_bar;
  double omega_l_tmp;
  arma::mat alpha_l_tmp;
  arma::mat beta_l_tmp;
  arma::vec eta_pl_tmp;
  double xi_star_bar;
  arma::mat omega_nl_tmp;
  arma::mat alpha_nl_tmp;
  arma::mat beta_nl_tmp;
  arma::rowvec xi_star_bar_nl;
  double sFct;
  double alpha_0_bar;
  double xi_star;
  arma::vec xi_starnl;
  ////////////////////////////////////////////////////
  //////////////////// Start MCMC ////////////////////
  ///////////////////////////////////////////////////
  
  for (int t = 0; t<iter; t++) {
    // update pi start linear
    pi_star_l = update_piSC(hyperpar(5), hyperpar(6), gamma_star_l, hyperpar(4));
    // update pi start non linear
    pi_star_nl = update_piSC(hyperpar(7), hyperpar(8), gamma_star_nl, hyperpar(4));
    // Time
    auto start_emLinear = std::chrono::high_resolution_clock::now();
    //////////////////// effect modifiers linear peNMIG ////////////////////
    for (int j = 0; j<(p-1); j++) {
      for (int k = (j+1); k<p; k++) {
        // if true I can update the interaction effect
        // weak heredity linear
        // bool wh = (gamma_0_l(j) != hyperpar(4)) || (gamma_0_l(k) != hyperpar(4));
        // if true I can update the interaction effect
        // voglio che per adesso valgano entrambi quindi se
        // almeno uno dei due è true allora calcolo gli effetti di interazione
        if (ha == 1) {
          if ((gamma_0_l(j) != hyperpar(4)) || (gamma_0_l(k) != hyperpar(4))) {
            // update gamma inclusion parameters
            gamma_star_l(j,k) = update_gammaScaC(pi_star_l, hyperpar(4), alpha_star_l(j,k), tau_star_l(j,k));
            // update tau variance
            tau_star_l(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_l(j,k), gamma_star_l(j,k));
            // update alpha star linear and non linear with a MH
            // proposed alpha
            alpha_star_bar = alpha_star_l(j,k) + R::rnorm(0, mht(0));
            // compute the new beta_l (temp element)
            omega_l_tmp = xi_star_l(j,k) * alpha_star_bar;
            arma::mat alpha_l_tmp = alpha_l;
            for (int jj = 0; jj<p; jj++) {
              alpha_l_tmp.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  if ((kn == k) && (jj == j)) {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l_tmp;
                  } else {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l(jj, kn);
                  }
                }
              }
            }
            beta_l_tmp = beta_l;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l_tmp(i, jj) = xi_l(jj) * alpha_l_tmp(i,jj);
              }
            }
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l_tmp, X_nl, beta_nl);
            // update linear alpha star
            // List
            List uasl = update_alphaC(y, sigma, tau_star_l(j,k), gamma_star_l(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_l(j,k));
            alpha_star_l(j,k) = uasl[0];
            int alpha_acc = uasl[1];
            alpha_star_l_acc(j,k) = alpha_star_l_acc(j,k) + alpha_acc;
            // compute linear omega
            omega_l(j,k) = alpha_star_l(j,k)*xi_star_l(j,k);
            // compute alpha linear
            for (int jj = 0; jj<p; jj++) {
              alpha_l.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  alpha_l.col(jj) = alpha_l.col(jj) + X_l.col(kn)*omega_l(jj, kn);
                }
              }
            }
            // compute linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l(i, jj) = xi_l(jj)*alpha_l(i, jj);
              }
            }
            // Linear Predictor with the proposed alpha
            eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
            // update m_star_l linear
            m_star_l(j,k) = update_mCsca(xi_star_l(j,k));
            // update xi star linear
            // proposed xi star
            xi_star_bar = xi_star_l(j,k) + R::rnorm(0, mht(1));
            // compute the new beta_l (temp element)
            omega_l_tmp =  alpha_star_l(j,k) * xi_star_bar;
            alpha_l_tmp = alpha_l;
            for (int jj = 0; jj<p; jj++) {
              alpha_l_tmp.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  if ((kn == k) && (jj == j)) {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l_tmp;
                  } else {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l(jj,kn);
                  }
                }
              }
            }
            beta_l_tmp = beta_l;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l_tmp(i,jj) = alpha_l_tmp(i,jj)*xi_l(jj);
              }
            }
            // Linear Predictor with the proposed xi star
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l_tmp, X_nl, beta_nl);
            // update xi linear star
            List uxsl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_star_l(j,k), xi_star_bar, xi_star_l(j,k));
            xi_star_l(j,k) = uxsl[0];
            int xiSacc = uxsl[1];
            xi_star_l_acc(j,k) = xi_star_l_acc(j,k) + xiSacc;
            omega_l(j,k) = alpha_star_l(j,k) * xi_star_l(j,k);
            // compute linear alpha
            for (int jj = 0; jj<p; jj++) {
              alpha_l.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  alpha_l.col(jj) = alpha_l.col(jj) + X_l.col(kn)*omega_l(jj, kn);
                }
              }
            }
            // compute linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l(i, jj) = alpha_l(i, jj)*xi_l(jj);
              }
            }
            // Linear Predictor with the proposed xi star
            eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
          }
        }
        
        if (ha == 2) {
          if ((gamma_0_l(j) != hyperpar(4)) && (gamma_0_l(k) != hyperpar(4))) {
            // update gamma inclusion parameters
            gamma_star_l(j,k) = update_gammaScaC(pi_star_l, hyperpar(4), alpha_star_l(j,k), tau_star_l(j,k));
            // update tau variance
            tau_star_l(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_l(j,k), gamma_star_l(j,k));
            // update alpha star linear and non linear with a MH
            // proposed alpha
            alpha_star_bar = alpha_star_l(j,k) + R::rnorm(0, mht(0));
            // compute the new beta_l (temp element)
            omega_l_tmp = xi_star_l(j,k) * alpha_star_bar;
            arma::mat alpha_l_tmp = alpha_l;
            for (int jj = 0; jj<p; jj++) {
              alpha_l_tmp.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  if ((kn == k) && (jj == j)) {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l_tmp;
                  } else {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l(jj, kn);
                  }
                }
              }
            }
            beta_l_tmp = beta_l;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l_tmp(i, jj) = xi_l(jj) * alpha_l_tmp(i,jj);
              }
            }
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l_tmp, X_nl, beta_nl);
            // update linear alpha star
            // List
            List uasl = update_alphaC(y, sigma, tau_star_l(j,k), gamma_star_l(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_l(j,k));
            alpha_star_l(j,k) = uasl[0];
            int alpha_acc = uasl[1];
            alpha_star_l_acc(j,k) = alpha_star_l_acc(j,k) + alpha_acc;
            // compute linear omega
            omega_l(j,k) = alpha_star_l(j,k)*xi_star_l(j,k);
            // compute alpha linear
            for (int jj = 0; jj<p; jj++) {
              alpha_l.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  alpha_l.col(jj) = alpha_l.col(jj) + X_l.col(kn)*omega_l(jj, kn);
                }
              }
            }
            // compute linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l(i, jj) = xi_l(jj)*alpha_l(i, jj);
              }
            }
            // Linear Predictor with the proposed alpha
            eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
            // update m_star_l linear
            m_star_l(j,k) = update_mCsca(xi_star_l(j,k));
            // update xi star linear
            // proposed xi star
            xi_star_bar = xi_star_l(j,k) + R::rnorm(0, mht(1));
            // compute the new beta_l (temp element)
            omega_l_tmp =  alpha_star_l(j,k) * xi_star_bar;
            alpha_l_tmp = alpha_l;
            for (int jj = 0; jj<p; jj++) {
              alpha_l_tmp.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  if ((kn == k) && (jj == j)) {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l_tmp;
                  } else {
                    alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l(jj,kn);
                  }
                }
              }
            }
            beta_l_tmp = beta_l;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l_tmp(i,jj) = alpha_l_tmp(i,jj)*xi_l(jj);
              }
            }
            // Linear Predictor with the proposed xi star
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l_tmp, X_nl, beta_nl);
            // update xi linear star
            List uxsl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_star_l(j,k), xi_star_bar, xi_star_l(j,k));
            xi_star_l(j,k) = uxsl[0];
            int xiSacc = uxsl[1];
            xi_star_l_acc(j,k) = xi_star_l_acc(j,k) + xiSacc;
            omega_l(j,k) = alpha_star_l(j,k) * xi_star_l(j,k);
            // compute linear alpha
            for (int jj = 0; jj<p; jj++) {
              alpha_l.col(jj) = alpha_0_l(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  alpha_l.col(jj) = alpha_l.col(jj) + X_l.col(kn)*omega_l(jj, kn);
                }
              }
            }
            // compute linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_l(i, jj) = alpha_l(i, jj)*xi_l(jj);
              }
            }
            // Linear Predictor with the proposed xi star
            eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
          }
        }
        
        if (ha == 0) {
          // update gamma inclusion parameters
          gamma_star_l(j,k) = update_gammaScaC(pi_star_l, hyperpar(4), alpha_star_l(j,k), tau_star_l(j,k));
          // update tau variance
          tau_star_l(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_l(j,k), gamma_star_l(j,k));
          // update alpha star linear and non linear with a MH
          // proposed alpha
          alpha_star_bar = alpha_star_l(j,k) + R::rnorm(0, mht(0));
          // compute the new beta_l (temp element)
          omega_l_tmp = xi_star_l(j,k) * alpha_star_bar;
          arma::mat alpha_l_tmp = alpha_l;
          for (int jj = 0; jj<p; jj++) {
            alpha_l_tmp.col(jj) = alpha_0_l(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                if ((kn == k) && (jj == j)) {
                  alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l_tmp;
                } else {
                  alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l(jj, kn);
                }
              }
            }
          }
          beta_l_tmp = beta_l;
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_l_tmp(i, jj) = xi_l(jj) * alpha_l_tmp(i,jj);
            }
          }
          // new linear predictor WITH the proposed alpha star
          eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l_tmp, X_nl, beta_nl);
          // update linear alpha star
          // List
          List uasl = update_alphaC(y, sigma, tau_star_l(j,k), gamma_star_l(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_l(j,k));
          alpha_star_l(j,k) = uasl[0];
          int alpha_acc = uasl[1];
          alpha_star_l_acc(j,k) = alpha_star_l_acc(j,k) + alpha_acc;
          // compute linear omega
          omega_l(j,k) = alpha_star_l(j,k)*xi_star_l(j,k);
          // compute alpha linear
          for (int jj = 0; jj<p; jj++) {
            alpha_l.col(jj) = alpha_0_l(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                alpha_l.col(jj) = alpha_l.col(jj) + X_l.col(kn)*omega_l(jj, kn);
              }
            }
          }
          // compute linear beta
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_l(i, jj) = xi_l(jj)*alpha_l(i, jj);
            }
          }
          // Linear Predictor with the proposed alpha
          eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
          // update m_star_l linear
          m_star_l(j,k) = update_mCsca(xi_star_l(j,k));
          // update xi star linear
          // proposed xi star
          xi_star_bar = xi_star_l(j,k) + R::rnorm(0, mht(1));
          // compute the new beta_l (temp element)
          omega_l_tmp =  alpha_star_l(j,k) * xi_star_bar;
          alpha_l_tmp = alpha_l;
          for (int jj = 0; jj<p; jj++) {
            alpha_l_tmp.col(jj) = alpha_0_l(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                if ((kn == k) && (jj == j)) {
                  alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l_tmp;
                } else {
                  alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(kn)*omega_l(jj,kn);
                }
              }
            }
          }
          beta_l_tmp = beta_l;
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_l_tmp(i,jj) = alpha_l_tmp(i,jj)*xi_l(jj);
            }
          }
          // Linear Predictor with the proposed xi star
          eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l_tmp, X_nl, beta_nl);
          // update xi linear star
          List uxsl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_star_l(j,k), xi_star_bar, xi_star_l(j,k));
          xi_star_l(j,k) = uxsl[0];
          int xiSacc = uxsl[1];
          xi_star_l_acc(j,k) = xi_star_l_acc(j,k) + xiSacc;
          omega_l(j,k) = alpha_star_l(j,k) * xi_star_l(j,k);
          // compute linear alpha
          for (int jj = 0; jj<p; jj++) {
            alpha_l.col(jj) = alpha_0_l(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                alpha_l.col(jj) = alpha_l.col(jj) + X_l.col(kn)*omega_l(jj, kn);
              }
            }
          }
          // compute linear beta
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_l(i, jj) = alpha_l(i, jj)*xi_l(jj);
            }
          }
          // Linear Predictor with the proposed xi star
          eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
        }
      } // end linear k
    } // end linear j

    auto stop_emLinear = std::chrono::high_resolution_clock::now();
    auto duration_emLinear = std::chrono::duration_cast<std::chrono::microseconds>(stop_emLinear - start_emLinear).count();

    // Time
    auto start_emNLinear = std::chrono::high_resolution_clock::now();
    //////////////////// effect modifiers non-linear peNMIG ////////////////////
    for (int j = 0; j<(p-1); j++) {
      for (int k = (j+1); k<p; k++) {
        
        if (ha == 1) {
          if ((gamma_0_l(j) != hyperpar(4)) || (gamma_0_l(k) != hyperpar(4))) {
            // update gamma star non linear
            gamma_star_nl(j,k) = update_gammaScaC(pi_star_nl, hyperpar(4), alpha_star_nl(j,k), tau_star_nl(j,k));
            // update tau star non linear
            tau_star_nl(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_nl(j,k), gamma_star_nl(j,k));
            // update alpha non linear star
            // proposed alpha
            alpha_star_bar = alpha_star_nl(j,k) + R::rnorm(0, mht(2));
            omega_nl_tmp = xi_star_nl(j, span(cd[k], cd[k+1]-1)) * alpha_star_bar;
            alpha_nl_tmp = alpha_nl;
            for (int jj = 0; jj<p; jj++) {
              alpha_nl_tmp.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  if ((kn == k) && (jj == j)) {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl_tmp.t();
                  } else {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                  }
                }
              }
            }
            // compute the beta non linear temp
            beta_nl_tmp = beta_nl;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl_tmp(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl_tmp(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl_tmp);
            // update alpha star non linear
            List uasnl = update_alphaC(y, sigma, tau_star_nl(j,k), gamma_star_nl(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_nl(j,k));
            alpha_star_nl(j,k) = uasnl[0];
            int acc_anl = uasnl[1];
            alpha_star_nl_acc(j,k) = alpha_star_nl_acc(j,k) + acc_anl;
            omega_nl(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1))*alpha_star_nl(j,k);
            // compute the linear predictor with the new alpha star non linear
            // compute the alpha non linear
            for (int jj = 0; jj<p; jj++) {
              alpha_nl.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  alpha_nl.col(jj) = alpha_nl.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                }
              }
            }
            // compute non linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // compute the linear predictor this the proposed alpha
            eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
            // update m star non linear vector
            m_star_nl(j, span(cd[k], cd[k+1]-1)) = update_mCvec(xi_star_nl(j, span(cd[k], cd[k+1]-1))).t();
            // update xi star non linear
            // proposed xi star non linear
            xi_star_bar_nl = xi_star_nl(j, span(cd[k], cd[k+1]-1)) + R::rnorm(0, mht(3));
            // compute the new beta 
            omega_nl_tmp = xi_star_bar_nl * alpha_star_nl(j,k);
            alpha_nl_tmp = alpha_nl;
            for (int jj = 0; jj<p; jj++) {
              alpha_nl_tmp.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  if ((kn == k) && (jj == j)) {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl_tmp.t();
                  } else {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                  }
                }
              }
            }
            // beta non linear temp
            beta_nl_tmp = beta_nl;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl_tmp(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl_tmp(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // compute the linear predictor this the proposed xi
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl_tmp);
            // update xi star non-linear
            List uxsnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_star_nl(j, span(cd[k], cd[k+1]-1)).t(), xi_star_bar_nl.t(), xi_star_nl(j, span(cd[k], cd[k+1]-1)).t());
            arma::rowvec resultXiSnl = uxsnl[0];
            xi_star_nl(j, span(cd[k], cd[k+1]-1)) = resultXiSnl;
            arma::rowvec accxisnl = uxsnl[1];
            xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) = xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) + accxisnl;
            omega_nl(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1))*alpha_star_nl(j,k);
            
            // compute the linear predictor with the proposed xi
            // compute the alpha non linear
            for (int jj = 0; jj<p; jj++) {
              alpha_nl.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  alpha_nl.col(jj) = alpha_nl.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t(); 
                }
              }
            }
            // compute the non linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // compute the linear predictor
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
          }
        }
        
        if (ha == 2) {
          if ((gamma_0_l(j) != hyperpar(4)) && (gamma_0_l(k) != hyperpar(4))) {
            // update gamma star non linear
            gamma_star_nl(j,k) = update_gammaScaC(pi_star_nl, hyperpar(4), alpha_star_nl(j,k), tau_star_nl(j,k));
            // update tau star non linear
            tau_star_nl(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_nl(j,k), gamma_star_nl(j,k));
            // update alpha non linear star
            // proposed alpha
            alpha_star_bar = alpha_star_nl(j,k) + R::rnorm(0, mht(2));
            omega_nl_tmp = xi_star_nl(j, span(cd[k], cd[k+1]-1)) * alpha_star_bar;
            alpha_nl_tmp = alpha_nl;
            for (int jj = 0; jj<p; jj++) {
              alpha_nl_tmp.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  if ((kn == k) && (jj == j)) {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl_tmp.t();
                  } else {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                  }
                }
              }
            }
            // compute the beta non linear temp
            beta_nl_tmp = beta_nl;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl_tmp(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl_tmp(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl_tmp);
            // update alpha star non linear
            List uasnl = update_alphaC(y, sigma, tau_star_nl(j,k), gamma_star_nl(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_nl(j,k));
            alpha_star_nl(j,k) = uasnl[0];
            int acc_anl = uasnl[1];
            alpha_star_nl_acc(j,k) = alpha_star_nl_acc(j,k) + acc_anl;
            omega_nl(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1))*alpha_star_nl(j,k);
            // compute the linear predictor with the new alpha star non linear
            // compute the alpha non linear
            for (int jj = 0; jj<p; jj++) {
              alpha_nl.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  alpha_nl.col(jj) = alpha_nl.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                }
              }
            }
            // compute non linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // compute the linear predictor this the proposed alpha
            eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
            // update m star non linear vector
            m_star_nl(j, span(cd[k], cd[k+1]-1)) = update_mCvec(xi_star_nl(j, span(cd[k], cd[k+1]-1))).t();
            // update xi star non linear
            // proposed xi star non linear
            xi_star_bar_nl = xi_star_nl(j, span(cd[k], cd[k+1]-1)) + R::rnorm(0, mht(3));
            // compute the new beta 
            omega_nl_tmp = xi_star_bar_nl * alpha_star_nl(j,k);
            alpha_nl_tmp = alpha_nl;
            for (int jj = 0; jj<p; jj++) {
              alpha_nl_tmp.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  if ((kn == k) && (jj == j)) {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl_tmp.t();
                  } else {
                    alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                  }
                }
              }
            }
            // beta non linear temp
            beta_nl_tmp = beta_nl;
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl_tmp(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl_tmp(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // compute the linear predictor this the proposed xi
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl_tmp);
            // update xi star non-linear
            List uxsnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_star_nl(j, span(cd[k], cd[k+1]-1)).t(), xi_star_bar_nl.t(), xi_star_nl(j, span(cd[k], cd[k+1]-1)).t());
            arma::rowvec resultXiSnl = uxsnl[0];
            xi_star_nl(j, span(cd[k], cd[k+1]-1)) = resultXiSnl;
            arma::rowvec accxisnl = uxsnl[1];
            xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) = xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) + accxisnl;
            omega_nl(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1))*alpha_star_nl(j,k);
            
            // compute the linear predictor with the proposed xi
            // compute the alpha non linear
            for (int jj = 0; jj<p; jj++) {
              alpha_nl.col(jj) = alpha_0_nl(jj)*vecOnes;
              if (jj != (p-1)) {
                for (int kn = (jj+1); kn<p; kn++) {
                  arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                  arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                  alpha_nl.col(jj) = alpha_nl.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t(); 
                }
              }
            }
            // compute the non linear beta
            for (int jj = 0; jj<p; jj++) {
              for (int i = 0; i<nobs; i++) {
                beta_nl(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
              }
            }
            // compute the linear predictor
            eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
          }
        }
        
        if (ha == 0) {
          // update gamma star non linear
          gamma_star_nl(j,k) = update_gammaScaC(pi_star_nl, hyperpar(4), alpha_star_nl(j,k), tau_star_nl(j,k));
          // update tau star non linear
          tau_star_nl(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_nl(j,k), gamma_star_nl(j,k));
          // update alpha non linear star
          // proposed alpha
          alpha_star_bar = alpha_star_nl(j,k) + R::rnorm(0, mht(2));
          omega_nl_tmp = xi_star_nl(j, span(cd[k], cd[k+1]-1)) * alpha_star_bar;
          alpha_nl_tmp = alpha_nl;
          for (int jj = 0; jj<p; jj++) {
            alpha_nl_tmp.col(jj) = alpha_0_nl(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                if ((kn == k) && (jj == j)) {
                  alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl_tmp.t();
                } else {
                  alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                }
              }
            }
          }
          // compute the beta non linear temp
          beta_nl_tmp = beta_nl;
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_nl_tmp(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl_tmp(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
            }
          }
          // new linear predictor WITH the proposed alpha star
          eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl_tmp);
          // update alpha star non linear
          List uasnl = update_alphaC(y, sigma, tau_star_nl(j,k), gamma_star_nl(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_nl(j,k));
          alpha_star_nl(j,k) = uasnl[0];
          int acc_anl = uasnl[1];
          alpha_star_nl_acc(j,k) = alpha_star_nl_acc(j,k) + acc_anl;
          omega_nl(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1))*alpha_star_nl(j,k);
          // compute the linear predictor with the new alpha star non linear
          // compute the alpha non linear
          for (int jj = 0; jj<p; jj++) {
            alpha_nl.col(jj) = alpha_0_nl(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                alpha_nl.col(jj) = alpha_nl.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
              }
            }
          }
          // compute non linear beta
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_nl(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
            }
          }
          // compute the linear predictor this the proposed alpha
          eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
          // update m star non linear vector
          m_star_nl(j, span(cd[k], cd[k+1]-1)) = update_mCvec(xi_star_nl(j, span(cd[k], cd[k+1]-1))).t();
          // update xi star non linear
          // proposed xi star non linear
          xi_star_bar_nl = xi_star_nl(j, span(cd[k], cd[k+1]-1)) + R::rnorm(0, mht(3));
          // compute the new beta 
          omega_nl_tmp = xi_star_bar_nl * alpha_star_nl(j,k);
          alpha_nl_tmp = alpha_nl;
          for (int jj = 0; jj<p; jj++) {
            alpha_nl_tmp.col(jj) = alpha_0_nl(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                if ((kn == k) && (jj == j)) {
                  alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl_tmp.t();
                } else {
                  alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t();
                }
              }
            }
          }
          // beta non linear temp
          beta_nl_tmp = beta_nl;
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_nl_tmp(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl_tmp(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
            }
          }
          // compute the linear predictor this the proposed xi
          eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl_tmp);
          // update xi star non-linear
          List uxsnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_star_nl(j, span(cd[k], cd[k+1]-1)).t(), xi_star_bar_nl.t(), xi_star_nl(j, span(cd[k], cd[k+1]-1)).t());
          arma::rowvec resultXiSnl = uxsnl[0];
          xi_star_nl(j, span(cd[k], cd[k+1]-1)) = resultXiSnl;
          arma::rowvec accxisnl = uxsnl[1];
          xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) = xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) + accxisnl;
          omega_nl(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1))*alpha_star_nl(j,k);
          
          // compute the linear predictor with the proposed xi
          // compute the alpha non linear
          for (int jj = 0; jj<p; jj++) {
            alpha_nl.col(jj) = alpha_0_nl(jj)*vecOnes;
            if (jj != (p-1)) {
              for (int kn = (jj+1); kn<p; kn++) {
                arma::vec col_ind = myRange(cd[kn], cd[kn+1]-1);
                arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
                alpha_nl.col(jj) = alpha_nl.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[kn], cd[kn+1]-1)).t(); 
              }
            }
          }
          // compute the non linear beta
          for (int jj = 0; jj<p; jj++) {
            for (int i = 0; i<nobs; i++) {
              beta_nl(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
            }
          }
          // compute the linear predictor
          eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
        }

      } // end non linear k
    } // end non linear j

    auto stop_emNLinear = std::chrono::high_resolution_clock::now();
    auto duration_emNLinear = std::chrono::duration_cast<std::chrono::microseconds>(stop_emNLinear - start_emNLinear).count();
    
    // rescaling
    for (int j = 0; j<(p-1); j++) {
      for (int k = (j+1); k<p; k++) {
        // linear
        sFct = arma::accu(abs(xi_star_l(j,k)));
        xi_star_l(j,k) = xi_star_l(j,k)/sFct;
        alpha_star_l(j,k) = sFct*alpha_star_l(j,k);
        // non linear
        sFct = d(k)/arma::accu(abs(xi_star_nl(j, span(cd[k], cd[k+1]-1))));
        xi_star_nl(j, span(cd[k], cd[k+1]-1)) = sFct*xi_star_nl(j, span(cd[k], cd[k+1]-1));
        alpha_star_nl(j,k) = alpha_star_nl(j,k)/sFct;
        // compute linear omega
        omega_l(j,k) = alpha_star_l(j,k)*xi_star_l(j,k);
        // compute non linear omega
        omega_nl(j, span(cd[k], cd[k+1]-1)) = alpha_star_nl(j,k)*xi_star_nl(j, span(cd[k], cd[k+1]-1));
      }
    }
    // update pi 0 linear
    pi_0_l = update_piNSC(hyperpar(5), hyperpar(6), gamma_0_l, hyperpar(4));
    // update pi 0 non linear
    pi_0_nl = update_piNSC(hyperpar(7), hyperpar(8), gamma_0_nl, hyperpar(4));
    // update gamma 0 linear
    gamma_0_l = update_gammaVecC(pi_0_l, hyperpar(4), alpha_0_l, tau_0_l);
    // update gamma 0 non linear
    gamma_0_nl = update_gammaVecC(pi_0_nl, hyperpar(4), alpha_0_nl, tau_0_nl);
    // update tau 0 linear
    for (int j = 0; j<p; j++) {
      tau_0_l(j) = update_tauC(hyperpar(0), hyperpar(1), alpha_0_l(j), gamma_0_l(j));
    }
    // update tau 0 non linear
    for (int j = 0; j<p; j++) {
      tau_0_nl(j) = update_tauC(hyperpar(0), hyperpar(1), alpha_0_nl(j), gamma_0_nl(j));
    }
    // update alpha 0 linear
    for (int j = 0; j<p; j++) {
      // proposed alpha
      alpha_0_bar = alpha_0_l(j) + R::rnorm(0, mht(4));
      // compute the new beta linear tmp
      alpha_l_tmp = alpha_l;
      for (int jj = 0; jj<p; jj++) {
        if (jj == j) {
          alpha_l_tmp.col(jj) = alpha_0_bar*vecOnes;
        } else {
          alpha_l_tmp.col(jj) = alpha_0_l(jj)*vecOnes;
        }
      }
      for (int jj = 0; jj<p; jj++) {
        if (jj != (p-1)) {
          for (int k = (jj+1); k<p; k++) {
            alpha_l_tmp.col(jj) = alpha_l_tmp.col(jj) + X_l.col(k)*omega_l(jj,k);
          }
        }
      }
      beta_l_tmp = beta_l;
      for (int jj = 0; jj<p; jj++) {
        for (int i = 0; i<nobs; i++) {
          beta_l_tmp(i,jj) = alpha_l_tmp(i,jj)*xi_l(jj);
        }
      }
      // new linear predictor with the proposed alpha 0 linear
      arma::vec eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l_tmp, X_nl, beta_nl);
      // update alpha 0 linear
      List ua0l = update_alphaC(y, sigma, tau_0_l(j), gamma_0_l(j), eta_pl_tmp, eta_pl, alpha_0_bar, alpha_0_l(j));
      alpha_0_l(j) = ua0l[0];
      int accAlpha0 = ua0l[1];
      alpha_0_l_acc(j) = alpha_0_l_acc(j) + accAlpha0;
      // new beta linear
      beta_l.col(j) = alpha_l.col(j)*xi_l(j);
    }
    // update alpha 0 non linear
    for (int j = 0; j<p; j++) {
      // proposed alpha
      alpha_0_bar = alpha_0_nl(j) + R::rnorm(0, mht(4));
      // compute the new beta linear tmp
      alpha_nl_tmp = alpha_nl;
      for (int jj = 0; jj<p; jj++) {
        if (jj == j) {
          alpha_nl_tmp.col(jj) = alpha_0_bar*vecOnes;
        } else {
          alpha_nl_tmp.col(jj) = alpha_0_nl(jj)*vecOnes;
        }
      }
      for (int jj = 0; jj<p; jj++) {
        if (jj != (p-1)) {
          for (int k = (jj+1); k<p; k++) {
            arma::vec col_ind = myRange(cd[k], cd[k+1]-1);
            arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
            alpha_nl_tmp.col(jj) = alpha_nl_tmp.col(jj) + X_nl.cols(ucol_ind)*omega_nl(jj, span(cd[k], cd[k+1]-1)).t();
          }
        }
      }
      beta_nl_tmp = beta_nl;
      for (int jj = 0; jj<p; jj++) {
        for (int i = 0; i<nobs; i++) {
          beta_nl_tmp(i, span(cd[jj], cd[jj+1]-1)) = alpha_nl_tmp(i,jj)*xi_nl(span(cd[jj], cd[jj+1]-1)).t();
        }
      }
      // new linear predictor with the proposed alpha 0 non linear
      eta_pl_tmp = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl_tmp);
      // update alpha 0 linear
      List ua0nl = update_alphaC(y, sigma, tau_0_nl(j), gamma_0_nl(j), eta_pl_tmp, eta_pl, alpha_0_bar, alpha_0_nl(j));
      alpha_0_nl(j) = ua0nl[0];
      int accAnlpha0 = ua0nl[1];
      alpha_0_nl_acc(j) = alpha_0_nl_acc(j) + accAnlpha0;
      // new beta linear
      beta_nl.col(j) = alpha_nl.col(j)*xi_nl(j);
    }
    // compute alpha linear
    for (int j = 0; j<p; j++) {
      alpha_l.col(j) = alpha_0_l(j)*vecOnes;
      if (j != (p-1)) {
        for (int k = (j+1); k<p; k++) {
          alpha_l.col(j) = alpha_l.col(j) + X_l.col(k)*omega_l(j,k);
        }
      }
    }
    //// compute alpha non linear
    for (int j = 0; j<p; j++) {
      alpha_nl.col(j) = alpha_0_nl(j)*vecOnes;
      if (j != (p-1)) {
        for (int k = (j+1); k<p; k++) {
          arma::vec col_ind = myRange(cd[k], cd[k+1]-1);
          arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
          alpha_nl.col(j) = alpha_nl.col(j) + X_nl.cols(ucol_ind)*omega_nl(j, span(cd[k], cd[k+1]-1)).t();
        }
      }
    }
    // update m linear
    for (int j = 0; j<p; j++) {
      m_l(j) = update_mCsca(xi_l(j));
    }
    // update m non linear
    for (int l = 0; l<q; l++) {
      m_nl(l) = update_mCsca(xi_nl(l));
    }
    // update xi linear
    for (int j = 0; j<p; j++) {
      // proposed xi
      xi_star = xi_l(j) + R::rnorm(0, mht(6));
      // compute the new beta linear
      beta_l_tmp = beta_l;
      beta_l_tmp.col(j) = alpha_l.col(j)*xi_star;
      // compute the new linear predictor with the proposed xi
      // % for the element wise multiplication
      arma::vec eta_pl_tmp = eta_pl - X_l.col(j)%beta_l.col(j) + X_l.col(j)%beta_l_tmp.col(j);
      // update xi linear
      List uxl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_l(j), xi_star, xi_l(j));
      xi_l(j) = uxl[0];
      int acc_xil = uxl[1];
      xi_l_acc(j) = xi_l_acc(j) + acc_xil;
    }
    // update xi non linear
    xi_starnl = xi_nl + as<arma::vec>(wrap(Rcpp::rnorm(q, 0, mht(7))));
    beta_nl_tmp = beta_nl;
    for (int j = 0; j<p; j++) {
      for (int i = 0; i<nobs; i++) {
        beta_nl_tmp(i, span(cd[j], cd[j+1]-1)) = xi_starnl(span(cd[j], cd[j+1]-1)).t()*alpha_nl(i,j);
      }
      arma::vec col_ind = myRange(cd[j], cd[j+1]-1);
      arma::uvec ucol_ind = arma::conv_to<arma::uvec>::from(col_ind);
      arma::mat sumEta = X_nl.cols(ucol_ind)%beta_nl.cols(ucol_ind) + X_nl.cols(ucol_ind)%beta_nl_tmp.cols(ucol_ind);
      arma::vec eta_pl_tmp = eta_pl - arma::sum(sumEta, 1); // apply function
      List uxnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_nl(span(cd[j], cd[j+1]-1)), xi_starnl(span(cd[j], cd[j+1]-1)), xi_nl(span(cd[j], cd[j+1]-1)));
      arma::vec resXnl = uxnl[0];
      xi_nl(span(cd[j], cd[j+1]-1)) = resXnl;
      arma::vec accXinl = uxnl[1];
      xi_nl_acc(span(cd[j], cd[j+1]-1)) = xi_nl_acc(span(cd[j], cd[j+1]-1)) + accXinl;
    }
    // rescale alpha and xi linear
    for (int j = 0; j<p; j++) {
      // scaling factor
      sFct = arma::accu(abs(xi_l(j)));
      // rescale linear xi
      xi_l(j) = xi_l(j)/sFct;
      // rescale alpha linear
      alpha_l.col(j) = sFct*alpha_l.col(j);
    }
    // rescale alpha and xi non linear
    for (int j = 0; j<p; j++) {
      // scaling factor
      sFct = d(j)/arma::accu(abs(xi_nl(span(cd[j], cd[j+1]-1))));
      // rescale non linear xi
      xi_nl(span(cd[j], cd[j+1]-1)) = sFct*xi_nl(span(cd[j], cd[j+1]-1));
      // rescale non linear alpha
      alpha_nl.col(j) = alpha_nl.col(j)/sFct;
    }
    // update beta linear and non linear after rescaling alpha and xi
    // beta linear
    for (int j = 0; j<p; j++) {
      for (int i = 0; i<nobs; i++) {
        beta_l(i, j) = alpha_l(i, j)*xi_l(j);
      }
    }
    // beta non linear
    for (int j = 0; j<p; j++) {
      for (int i = 0; i<nobs; i++) {
        beta_nl(i, span(cd[j], cd[j+1]-1)) = alpha_nl(i,j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
      }
    }
    // compute the linear predictor
    arma::vec eta_pl = compLinPred(nobs, p, cd, eta0, X_l, beta_l, X_nl, beta_nl);
    // update intercept
    // linear predictor without intercept
    arma::vec eta_noInt = eta_pl - eta0*vecOnes;
    eta0 = updateInterceptC(y, nobs, eta_noInt, sigma);
    // update linear predictor
    eta_pl = eta_noInt + eta0*vecOnes;
    // update sigma variance
    sigma = update_sigmaC(y, eta_pl, hyperpar(2), hyperpar(3), nobs);
    // log-likelihood
    double logLik = arma::accu(dnormLogVec(y, eta_pl, sqrt(sigma)));

    // Time
    arma::vec timeVector(2);
    timeVector(0) = duration_emLinear;
    timeVector(1) = duration_emNLinear;
    // store resutls
    if(t%thin == 0 && t > burnin-1) { // we start from 0
      PI_S_l(idx) = pi_star_l;
      PI_S_nl(idx) = pi_star_nl;
      PI_0_l(idx) = pi_0_l;
      PI_0_nl(idx) = pi_0_nl;
      ETA_PL.row(idx) = eta_pl.t();
      ETA0(idx) = eta0;
      SIGMA(idx) = sigma;
      LOGLIKELIHOOD(idx) = logLik;
      GAMMA_S_l[idx] = gamma_star_l;
      GAMMA_S_nl[idx] = gamma_star_nl;
      TAU_S_l[idx] = tau_star_l;
      TAU_S_nl[idx] = tau_star_nl;
      ALPHA_S_l[idx] = alpha_star_l;
      ALPHA_S_nl[idx] = alpha_star_nl;
      M_S_l[idx] = m_star_l;
      M_S_nl[idx] = m_star_nl;
      XI_S_l[idx] = xi_star_l;
      XI_S_nl[idx] = xi_star_nl;
      OMEGA_l[idx] = omega_l;
      OMEGA_nl[idx] = omega_nl;
      M_l.row(idx) = m_l.t();
      M_nl.row(idx) = m_nl.t();
      XI_l.row(idx) = xi_l.t();
      XI_nl.row(idx) = xi_nl.t();
      BETA_l[idx] = beta_l;
      BETA_nl[idx] = beta_nl;
      GAMMA_0_l.row(idx) = gamma_0_l.t();
      GAMMA_0_nl.row(idx) = gamma_0_nl.t();
      TAU_0_l.row(idx) = tau_0_l.t();
      TAU_0_nl.row(idx) = tau_0_nl.t();
      ALPHA_0_l.row(idx) = alpha_0_l.t();
      ALPHA_0_nl.row(idx) = alpha_0_nl.t();
      ALPHA_l[idx] = alpha_l;
      ALPHA_nl[idx] = alpha_nl;
      // Time
      Prof_Mat.row(idx) = timeVector.t();
      idx = idx + 1;
    }
  } // end iteration
  ////////////////////////////////////////////////////
  //////////////////// End MCMC //////////////////////
  ///////////////////////////////////////////////////
  // Time 
  //auto stop = std::chrono::high_resolution_clock::now();
  //auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  
  // std::cout << "The computational time for the entire MCMC is " << duration/1000000;
  
  return List::create(Named("d") = d,
                      Named("intercept") = ETA0,
                      Named("linear_predictor") = ETA_PL,
                      Named("alpha_0_l") = ALPHA_0_l,
                      Named("alpha_0_nl") = ALPHA_0_nl,
                      Named("alpha_star_l") = ALPHA_S_l,
                      Named("alpha_star_nl") = ALPHA_S_nl,
                      //Named("m_l") = M_l,
                      //Named("m_nl") = M_nl,
                      // Named("m_star_l") = M_S_l,
                      // Named("m_star_nl") = M_S_nl,
                      Named("xi_l") = XI_l,
                      Named("xi_nl") = XI_nl,
                      Named("xi_star_l") = XI_S_l,
                      Named("xi_star_nl") = XI_S_nl,
                      //Named("tau_0_l") = TAU_0_l,
                      //Named("tau_0_nl") = TAU_0_nl,
                      //Named("tau_star_l") = TAU_S_l,
                      //Named("tau_star_nl") = TAU_S_nl
                      Named("gamma_0_l") = GAMMA_0_l,
                      Named("gamma_0_nl") = GAMMA_0_nl,
                      Named("gamma_star_l") = GAMMA_S_l,
                      Named("gamma_star_nl") = GAMMA_S_nl,
                      //Named("pi_0_l") = PI_0_l
                      //Named("pi_0_nl") = PI_0_nl,
                      //Named("pi_star_l") = PI_S_l,
                      Named("pi_star_nl") = PI_S_nl,
                      Named("sigma") = SIGMA,
                      Named("LogLikelihood") = LOGLIKELIHOOD,
                      //Named("acc_a_s_l") = alpha_star_l_acc/iter,
                      //Named("acc_a_s_nl") = alpha_star_nl_acc/iter,
                      //Named("acc_xi_s_l") = xi_star_l_acc/iter,
                      //Named("acc_xi_s_nl") = xi_star_nl_acc/iter,
                      //Named("acc_a_0_l") = alpha_0_l_acc/iter, 
                      //Named("acc_a_0_nl") = alpha_0_nl_acc/iter, 
                      //Named("acc_xi_l") = xi_l_acc/iter, 
                      //Named("acc_xi_nl") = xi_nl_acc/iter
                      //Named("Execution Time") = duration/1000000
                      Named("Profiling") = Prof_Mat
  );
}