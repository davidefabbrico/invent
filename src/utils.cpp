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

// Function to generate 'n' random numbers from a Gamma(shape, scale) distribution.
vec callrgamma(int n, double shape, double scale) {
  return vec(rgamma(n, shape, scale));
}

// Function to generate a random number from a normal distribution with  
// a given mean and standard deviation.  
// Uses randn(), which samples from a standard normal distribution (mean = 0, stddev = 1),  
// then scales and shifts the value to match the specified mean and standard deviation. 
double generate_normal(double mean, double stddev) {
  return mean + stddev * randn();
}

// INTERCEPT
double updateInterceptC(vec y, int nobs, vec lp_noInt, double sigma) {
  // Function to update the intercept term. Returns the updated intercept value.  
  double intercept = generate_normal(accu(y-lp_noInt)/(nobs+1), sigma/(nobs+1));
  return intercept;
}

// M scalar
double update_mCsca(double xi) {
  // Function to update the m vector. Returns the updated m value as a scalar.
  double m = 1.0/(1.0+exp(-2.0*xi));
  return m;
}

// M vector
vec update_mCvec(rowvec xi) {
  // Function to update the m vector. Returns the updated m value as a vector.
  int n = size(xi)(1);
  vec m(n);
  for (int i = 0; i<n; i++) {
    m(i) = 1.0/(1.0+exp(-2.0*xi(i)));
  } 
  return m;
} 

// TAU
double update_tauC(double at, double bt, double alpha, double gamma) {
  // Function to update the tau value. Returns the updated tau value.
  double tau = 1.0/(callrgamma(1, at+0.5, 1.0/(bt + pow(alpha, 2)/(2.0*gamma)))(0));
  return tau;
}

// DELTA F Non sctructured
// Function to count the occurrences of a specific value in a given Armadillo vector.  
// Iterates through the vector 'x' and increments the counter 'count' each time  
// an element matches the specified value 'val'. Returns the total count of occurrences.  
int deltaFNSC(vec x, double val) {
  int count = 0;
  int n = x.n_elem;
  for (int i = 0; i < n; i++) {
    if (x[i] == val) {
      count = count + 1;
    } 
  }
  return count;
} 

// DELTA F
// Function to count the occurrences of a specific value in the upper triangular  
// part (excluding the diagonal) of a given Armadillo matrix.  
// Iterates over all unique element pairs (j, k) where j < k and increments  
// the counter 'count' each time an element matches the specified value 'val'.  
// Returns the total count of occurrences.  
int deltaFSC(mat x, double val) {
  int count = 0;
  int n = size(x)(0);
  for (int j = 0; j<(n-1); j++) {
    for (int k = (j+1); k<n;k++) {
      if (x(j, k) == val) {
        count = count + 1;
      }
    } 
  }
  return count;
} 

// PI Non structured
double update_piNSC(double ap, double bp, vec gamma, double v0) {
  // Function to update the pi value (vector gamma). Returns the updated pi value.
  double pi = rbeta(1, ap + deltaFNSC(gamma, 1), bp + deltaFNSC(gamma, v0))(0);
  return pi;
} 

// PI Structured
double update_piSC(double ap, double bp, mat gamma, double v0) {
  // Function to update the pi value (matrix gamma). Returns the updated pi value.
  double pi = rbeta(1, ap + deltaFSC(gamma, 1), bp + deltaFSC(gamma, v0))(0);
  return pi;
} 

// SIGMA
double update_sigmaC(colvec y, colvec eta_pl, double as, double bs, int nobs) {
  // Function to update the sigma value. Returns the updated sigma value.
  double sigma = 1.0/(callrgamma(1, as + nobs/2.0, 1.0/(bs + accu(pow((y - eta_pl), 2))/2.0))(0));
  return sigma;
} 

// Function to update the value of gamma based on the given parameters pi, v0, alpha, and tau.  
// The function first computes two intermediate variables, d1 and d2, based on the input parameters.  
// It then calculates 'out', which is used to determine the probability (prob) of gamma being updated to 1.  
// If 'out' is very large or very small, the function sets gamma_out to 1 or v0, respectively.  
// Otherwise, it uses a random sample to decide whether gamma_out should be updated to 1 or remain at v0.  
// Returns the updated value of gamma_out as a scalar.  
double update_gammaScaC(double pi, double v0, double alpha, double tau) {
  // Update the GAMMA Scalar
  double gamma_out = 0;
  double prob = 0.0;
  double d1 = log(pi/(1.0-pi)) + log(v0) * 1.0/2.0;
  double d2 = (1.0 - v0) / (2.0 * v0);
  double out = exp(d1 + d2 * pow(alpha, 2)/tau);
  if (out > 10000) {
    prob = 1.0;
    gamma_out = 1;
  } else {
    if (out < 0.00001) {
      prob = 0.0;
      gamma_out = v0;
    } else { 
      prob = out / (1.0 + out);
      double samp = randu<vec>(1)(0); 
      if (samp < prob) {
        gamma_out = 1;
      } else { 
        gamma_out = v0;
      } 
    }
  }
  return gamma_out;
} 


// Update GAMMA Vector
// Function to update the value of gamma based on the given parameters pi, v0, alpha, and tau.  
// The function first computes two intermediate variables, d1 and d2, based on the input parameters.  
// It then calculates 'out', which is used to determine the probability (prob) of gamma being updated to 1.  
// If 'out' is very large or very small, the function sets gamma_out to 1 or v0, respectively.  
// Otherwise, it uses a random sample to decide whether gamma_out should be updated to 1 or remain at v0.  
// Returns the updated value of gamma_out as a vector.
vec update_gammaVecC(double pi, double v0, vec alpha, vec tau) {
  int p = alpha.n_elem;
  vec gamma_out(p);
  double prob = 0.0;
  double d1 = log(pi/(1.0-pi)) + log(v0) * 1.0/2.0;
  double d2 = (1.0 - v0) / (2.0 * v0);
  for (int j = 0; j<p; j++) {
    double out = exp(d1 + d2 * pow(alpha(j), 2)/tau(j));
    if (out > 10000) {
      prob = 1.0;
      gamma_out(j) = 1;
    } else {
      if (out < 0.00001) {
        prob = 0.0;
        gamma_out(j) = v0;
      } else { 
        prob = out / (1.0 + out);
        double samp = randu<vec>(1)(0);
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

// Function to compute the log of the probability density function (PDF) of a normal distribution.  
// The function takes the value 'x', the mean 'mean', and the standard deviation 'sd' as inputs.  
// It calculates the log of the normal distribution's PDF using the formula:
// log(PDF) = -0.5 * (log(2 * pi) + 2 * log(sd) + ((x - mean) / sd)^2)  
// Returns the computed log-PDF value for the given inputs.  
double normal_log_pdf(double x, double mean, double sd) {
  return -0.5 * (log(2 * M_PI) + 2 * log(sd) + pow((x - mean) / sd, 2));
}

// Function to compute the log of the probability density function (log-PDF) of a normal distribution  
// for each element in a given vector 'x', using corresponding means from the 'means' vector and a shared standard deviation 'sds'.  
// For each element x(i), the log-PDF is calculated using the 'normal_log_pdf' function.  
// The result is stored in a new vector 'res', which is returned as the output.  
vec dnormLogVec(vec x, vec means, double sds) {
  int n = x.n_elem;
  vec res(n);
  for(int i = 0; i < n; i++) {
    res(i) = normal_log_pdf(x(i), means(i), sds);
  }
  return res;
}

// ALPHA
// Function to update the value of the parameter alpha using a Metropolis-Hastings step.  
// The function computes the log-likelihood for the current and proposed values of alpha and uses them to calculate the acceptance ratio.  
// It then compares the ratio with a random sample from a uniform distribution to decide whether to accept or reject the new alpha value.  
// The function returns a List containing the updated value of alpha and the acceptance indicator (1 if accepted, 0 otherwise). 
List update_alphaC(vec y, double sigma, double tau, double gamma, vec eta_star, vec eta_pl, double alpha_star, double alpha) {
  int acc_alpha = 0;
  double a = accu(dnormLogVec(y, eta_star, sqrt(sigma)));
  double b = normal_log_pdf(alpha_star, 0, sqrt(tau*gamma));
  double num_alpha = a+b;
  double den_alpha = accu(dnormLogVec(y, eta_pl, sqrt(sigma))) + normal_log_pdf(alpha, 0, sqrt(tau*gamma));
  double ratio_alpha = num_alpha - den_alpha;
  double lsamp = log(randu<vec>(1)(0));
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
// Function to update the value of the parameter xi using a Metropolis-Hastings step.  
// The function computes the log-likelihood for the current and proposed values of xi and calculates the acceptance ratio.  
// It then compares the ratio with a random sample from a uniform distribution to decide whether to accept or reject the new xi value.  
// The function returns a List containing the updated value of xi and the acceptance indicator (1 if accepted, 0 otherwise).  
List update_xiLC(vec y, colvec eta_star, colvec eta_pl, double sigma, double m, double xi_star, double xi) {
  int acc_xi = 0;
  double num_xi = accu(dnormLogVec(y, eta_star, sqrt(sigma))) + normal_log_pdf(xi_star, m, 1);
  double den_xi = accu(dnormLogVec(y, eta_pl, sqrt(sigma))) + normal_log_pdf(xi, m, 1);
  double ratio_xi = num_xi - den_xi;
  double lsamp = log(randu<vec>(1)(0));
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
// Function to update the parameter xi (vector) using a Metropolis-Hastings step for multiple elements.  
// The function computes the log-likelihood for the current and proposed values of xi (as vectors) and calculates the acceptance ratio.  
// It then compares the ratio with a random sample from a uniform distribution to decide whether to accept or reject the new xi values.  
// The function returns a List containing the updated value of xi (as a vector) and the acceptance indicator for each element of xi (1 if accepted, 0 otherwise).  
List update_xiNLC(vec y, vec eta_star, vec eta_pl, double sigma, vec m, vec xi_star, vec xi) {
  int dj = xi_star.n_elem;
  vec acc_xi(dj);
  double num_xi = accu(dnormLogVec(y, eta_star, sqrt(sigma))) + accu(dnormLogVec(xi_star, m, 1));
  double den_xi = accu(dnormLogVec(y, eta_pl, sqrt(sigma))) + accu(dnormLogVec(xi, m, 1));
  double ratio_xi = num_xi - den_xi;
  double lsamp = log(randu<vec>(1)(0));
  vec xi_new;
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
// Function to compute the sign of a number. Returns 1 if x is positive, 
// -1 if x is negative. If x is 0, the function returns 0.
int mysign(double x) {
  if (x < 0) {
    return -1;
  } else if (x > 0) {
    return 1;
  } else { 
    return 0;
  }
} 


// Compute Linear Predictor
vec compLinPred(int nobs, double eta0, const mat& X_l, const mat& beta_l, 
                const mat& X_nl, const mat& beta_nl) {
  vec eta_pl(nobs, fill::value(eta0));
  vec sum_linear(nobs, fill::zeros);
  for (uword k = 0; k < X_l.n_cols; ++k) {
    sum_linear += X_l.col(k) % beta_l.col(k);
  }
  vec sum_nonlinear(nobs, fill::zeros);
  for (uword k = 0; k < X_nl.n_cols; ++k) {
    sum_nonlinear += X_nl.col(k) % beta_nl.col(k);
  }
  eta_pl += sum_linear + sum_nonlinear;
  return eta_pl;
}

// Initializes a p x p matrix with zeros and sets the values above the main diagonal to `hyperpar`
mat initGammaStar(int p, double hyperpar) {
  mat parameter(p, p, fill::zeros);
  for (uword j = 0; j < p; ++j) {
    const uword start_k = j + 1;
    if (start_k < p) {
      parameter(span(j), span(start_k, p-1)).fill(hyperpar); // Fills the upper part with hyperpar
    }
  }
  return parameter;
} 

// Initializes a p x p matrix with zeros and assigns values to the upper triangular part
// The values are sampled from a gamma distribution and then inverted (1/gamma)
mat initTauStar(int p, double a, double b) {
  mat parmeter(p, p, fill::zeros);
  int n = p*(p-1)/2; 
  vec gamma_samples = callrgamma(n, a, b); 
  parmeter.elem(trimatu_ind(size(parmeter), 1)) = 1.0 / gamma_samples;
  return parmeter;
}
 
// Initializes a vector of size p with values sampled from a gamma distribution
// The values are inverted (1/gamma)
vec initTauMain(int p, double a, double b) {
  vec gamma_samples = callrgamma(p, a, b); 
  return 1.0 / gamma_samples;
} 

// Initializes an nlp x q matrix with zeros and fills the upper part with 1.0
mat initXiStar(int nlp, int q) {
  mat parameter(nlp, q, fill::zeros);
  for (uword j = 0; j < min(nlp, q-1); ++j) {
    const int start_col = j + 1;
    if (start_col < q) {
      parameter.submat(span(j), span(start_col, q-1)).fill(1.0); // Fills with 1.0
    }
  } 
  return parameter;
} 

// Initializes a p x p matrix with zeros and assigns random values (+1 or -1) to the upper triangular part
mat initMStarL(int p) {
  mat parameter(p, p, fill::zeros); 
  mat rand_values = 2.0 * randu<mat>(p, p) - 1.0; // Matrix with random values between -1 and 1
  parameter(trimatu_ind(size(parameter), 1)) = sign(rand_values(trimatu_ind(size(parameter), 1)));
  return parameter;
} 

// Initializes a vector of length p with random values +1 or -1
vec initMMain(int p) {
  vec rand_values = 2.0 * randu<vec>(p) - 1.0; // Vector with random values between -1 and 1
  return sign(rand_values);
} 

// Initializes an nlp x q matrix with zeros and assigns values +1 or -1 to the upper triangular part
mat initMStarNL(int nlp, int q) {
  mat m_star_nl(nlp, q, fill::zeros); 
  mat rand_values = 2.0 * randu<mat>(nlp, q) - 1.0; // Matrix with random values between -1 and 1
  for (uword j = 0; j < nlp; ++j) {
    for (uword k = j + 1; k < q; ++k) {
      m_star_nl(j, k) = (rand_values(j, k) > 0) ? 1.0 : -1.0; // Assigns +1 or -1 based on the sign of the random value
    }
  } 
  return m_star_nl;
} 

// Body MCMC
// [[Rcpp::export]]
List bodyMCMC(vec y, int p, int nobs, vec cd, vec d, mat X_l, 
              mat X_nl, mat X_val_l, mat X_val_nl, vec hyperpar, 
              vec mht, int n_cat, int iter, int burnin, int thin, int ha, 
              bool detail = false, bool pb = true) {
  // Body MCMC
  // Function to perform a Markov Chain Monte Carlo (MCMC) process to sample from a posterior distribution.
  // The function uses the provided input data and hyperparameters to run the MCMC algorithm for the specified number of iterations, 
  // with options for burn-in and thinning to improve sampling efficiency. It also includes options for tracking progress and printing detailed results.
  //
  // **Inputs:**
  // - `y`: A vector containing the observed data (dependent variable) for the model.
  // - `p`: The number of covariates.
  // - `nobs`: The total number of observations (data points) in the dataset.
  // - `cd`: A vector representing the cumulative sum of the number of spline bases for the covariate.
  // - `d`: A vector indicating the number of spline bases for each covariate.
  // - `X_l`: A matrix of explanatory variables (covariates) for the linear component of the model.
  // - `X_nl`: A matrix of explanatory variables (covariates) for the non-linear component of the model.
  // - `X_val_l`: A matrix of validation covariates for the linear component.
  // - `X_val_nl`: A matrix of validation covariates for the non-linear component.
  // - `hyperpar`: A vector of hyperparameters that define the prior distribution for the model parameters.
  // - `mht`: A vector containing Metropolis-Hastings-related parameters for proposing new values during sampling.
  // - `n_cat`: The number of categorical variables in the dataset.
  // - `iter`: The total number of iterations for the MCMC sampling process.
  // - `burnin`: The number of initial iterations to discard during the burn-in period.
  // - `thin`: The thinning interval, specifying how often to retain a sample.
  // - `ha`: The heredity assumption for the model, potentially defining relationships between variables.
  // - `detail`: A boolean flag (default: false) indicating whether to retain all MCMC chain samples or not.
  // - `pb`: A boolean flag (default: true) specifying whether to display a progress bar during sampling.
  if (thin <= 0) {
    cout << "Look at the thin value, it must be greater than 0" << "\n";
  }
  // Time 
  auto start = chrono::high_resolution_clock::now();
  ////////////////////////////////////////////////////
  ////////////////// Initial value //////////////////
  ///////////////////////////////////////////////////
  // no assumption = 0
  // weak heredity = 1
  // strong heredity = 2
  int q = accu(d);
  // Number of non linear basis (n-n_cat)
  int nlp = p - n_cat; // the categorical one are the last ones
  vec vecOnes = ones(nobs, 1);
  // pi star linear 
  double pi_star_l = rbeta(1, 1, 1)(0);
  // pi star non linear
  double pi_star_nl = rbeta(1, 1, 1)(0);
  // pi 0 linear
  double pi_0_l = rbeta(1, 1, 1)(0);
  // pi 0 non linear
  double pi_0_nl = rbeta(1, 1, 1)(0);
  // save hyperpar4 to avoid accessing it repeatedly
  const double hyperpar4 = hyperpar(4);
  // gamma star linear
  // upper trinagular matrix
  mat gamma_star_l = initGammaStar(p, hyperpar4);
  // gamma star non linear
  // upper triangular matrix 
  mat gamma_star_nl = initGammaStar(nlp, hyperpar4);
  // gamma 0 linear (vector of dimension p)
  vec gamma_0_l(p, fill::value(hyperpar4));
  // gamma 0 non linear (vector of dimension nlp)
  vec gamma_0_nl(nlp, fill::value(hyperpar4));
  // tau star linear
  mat tau_star_l = initTauStar(p, 1.0, 1.0);
  // tau star non linear
  mat tau_star_nl = initTauStar(nlp, 1.0, 1.0);
  // tau linear (variance vector)
  vec tau_0_l = initTauMain(p, 1.0, 1.0);
  // tau non linear (variance vector)
  vec tau_0_nl = initTauMain(nlp, 1.0, 1.0);
  // alpha star linear square matrix of dimension p (for interaction)
  mat alpha_star_l(p, p, fill::zeros);
  // Non linear
  mat alpha_star_nl(nlp, nlp, fill::zeros);
  // xi star linear (square upper trinangular matrix of dimension p)
  mat xi_star_l = initGammaStar(p, 1.0);
  // xi star non linear, matrix of dimension pxq where q is the sum of the basis
  mat xi_star_nl = initXiStar(nlp, q);
  // omega linear
  mat omega_l(p,p);
  // omega non linear
  mat omega_nl(nlp,q);
  // alpha 0 linear vector (theta in thesis)
  vec alpha_0_l(p);
  // alpha 0 non linear vector (theta in thesis)
  vec alpha_0_nl(nlp);
  // alpha linear
  mat alpha_l(nobs, p);
  // alpha non linear
  mat alpha_nl(nobs, nlp);
  // m star linear matrix of dimension pxp
  mat m_star_l = initMStarL(p);
  // m star non linear 
  mat m_star_nl = initMStarNL(nlp, q);
  // m linear vector
  vec m_l = initMMain(p);
  // m non linear vector
  vec m_nl = initMMain(q);
  // xi linear
  vec xi_l = ones<vec>(p);
  // xi non linear
  vec xi_nl = ones<vec>(q);
  // beta linear
  mat beta_l(nobs, p);
  // beta non linear
  mat beta_nl(nobs, q);
  // intercept
  double eta0 = randn();
  // variance of the normal model
  double sigma = 1.0/callrgamma(1, 1, 1)(0);
  // Linear Predictor
  vec eta_pl = compLinPred(nobs, eta0, X_l, beta_l, X_nl, beta_nl);
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
  mat alpha_star_l_acc(p,p);
  mat alpha_star_nl_acc(nlp,nlp);
  // xi star acc
  mat xi_star_l_acc(p,p);
  mat xi_star_nl_acc(nlp,q);
  // alpha 0 (theta in thesis) acc
  vec alpha_0_l_acc(p);
  vec alpha_0_nl_acc(nlp);
  // xi acc
  vec xi_l_acc(p);
  vec xi_nl_acc(q);
  // pi star linear and non linear
  vec PI_S_l(nout);
  vec PI_S_nl(nout);
  // pi 0 linear and non linear
  vec PI_0_l(nout);
  vec PI_0_nl(nout);
  // linear predictor
  mat ETA_PL(nout, nobs);
  // intercept
  vec ETA0(nout);
  // sigma variance model
  vec SIGMA(nout);
  // log-likelihood
  vec LOGLIKELIHOOD(nout);
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
  mat M_l(nout, p);
  mat M_nl(nout, q);
  // xi linear and non linear
  mat XI_l(nout, p);
  mat XI_nl(nout, q);
  // Beta linear and non linear
  List BETA_l(nout);
  List BETA_nl(nout);
  // Gamma star linear and non linear
  mat GAMMA_0_l(nout, p);
  mat GAMMA_0_nl(nout, nlp);
  // Tau 0 linear and non linear
  mat TAU_0_l(nout, p);
  mat TAU_0_nl(nout, nlp);
  mat ALPHA_0_l(nout, p);
  mat ALPHA_0_nl(nout, nlp);
  // Alpha linear and non linear
  List ALPHA_l(nout);
  List ALPHA_nl(nout);
  // Init Parameters
  double alpha_star_bar;
  mat alpha_l_tmp = alpha_l;
  mat beta_l_tmp = beta_l;
  vec eta_pl_tmp;
  double xi_star_bar;
  mat omega_l_tmp = omega_l;
  mat omega_nl_tmp = omega_nl;
  mat alpha_nl_tmp = alpha_nl;
  mat beta_nl_tmp = beta_nl;
  rowvec xi_star_bar_nl;
  double sFct;
  double alpha_0_bar;
  double xi_star;
  vec xi_starnl;
  // Predictive
  int n_val = 0;
  if (!X_val_l.is_empty()) {
    n_val = X_val_l.n_rows;
  }
  mat Y_TILDE(nout, n_val);
  // predictive
  mat alpha_val_l(n_val, p);
  mat alpha_val_nl(n_val, nlp);
  mat beta_val_l(n_val, p);
  mat beta_val_nl(n_val, q);
  vec eta_pl_val(n_val);
  vec y_tilde(n_val);
  vec vecOnesVal = ones(n_val, 1);
  // detail = false, we choose to limit the growing of the interaction chain
  // in particuar the matrix.
  vec gamma_0_l_m(p); // main linear
  vec gamma_0_nl_m(nlp); // main nonlinear
  mat gamma_star_l_m(p,p); // interaction linear
  mat gamma_star_nl_m(nlp,nlp); // interaction nonlinear
  vec eta_pl_m(nobs); // linear predictor
  vec y_tilde_m(n_val); // prediction
  ////////////////////////////////////////////////////
  //////////////////// Start MCMC ////////////////////
  ///////////////////////////////////////////////////
  if (pb == true) {
    std::cout << "Running MCMC loop:\n";
  }
  const double hyperpar4_sq = hyperpar4 * hyperpar4;
  uvec active_j;
  umat active_pairs;
  for (int t = 0; t<iter; t++) {
    // update pi start linear
    pi_star_l = update_piSC(hyperpar(9), hyperpar(10), gamma_star_l, hyperpar4);
    // update pi start non linear
    pi_star_nl = update_piSC(hyperpar(11), hyperpar(12), gamma_star_nl, hyperpar4);
    // update gamma linear
    gamma_0_l = update_gammaVecC(pi_0_l, hyperpar4, alpha_0_l, tau_0_l);
    // check interaction selection indicators
    if (ha == 2) {
      active_j = find(gamma_0_l == hyperpar4);
      for (uword j : active_j) {
        if (j < p - 1) {
          gamma_star_l.submat(j, j+1, j, p-1).fill(hyperpar4);
        }
      }
    } else if (ha == 1) {
      active_pairs = conv_to<umat>::from(
        (gamma_0_l * gamma_0_l.t()) == hyperpar4_sq
      );
      active_pairs = trimatu(active_pairs, 1);
      gamma_star_l.elem(find(active_pairs > 0)).fill(hyperpar4);
    }
    //////////////////// effect modifiers linear peNMIG ////////////////////
    for (int j = 0; j<p; j++) {
      for (int k = (j+1); k<p; k++) {
        if (ha == 1) {
          if ((gamma_0_l(j) != hyperpar4) || (gamma_0_l(k) != hyperpar4)) {
            // update gamma inclusion parameters
            gamma_star_l(j,k) = update_gammaScaC(pi_star_l, hyperpar4, alpha_star_l(j,k), tau_star_l(j,k));
            // update tau variance
            tau_star_l(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_l(j,k), gamma_star_l(j,k));
            // update alpha star linear and non linear with a MH
            // proposed alpha
            alpha_star_bar = alpha_star_l(j,k) + mht(0)*randn(); // generate_normal(0, mht(0));
            // compute the new beta_l (temp element)
            omega_l_tmp = omega_l;
            omega_l_tmp(j,k) = xi_star_l(j,k) * alpha_star_bar;
            // change the j,k
            alpha_l_tmp.col(j) = alpha_0_l(j)*vecOnes;
            for (int kn = (j+1); kn<p; kn++) {
              alpha_l_tmp.col(j) = alpha_l_tmp.col(j) + X_l.col(kn)*omega_l_tmp(j,kn);
            }
            beta_l_tmp = beta_l;
            // compute linear beta
            beta_l_tmp.col(j) = alpha_l_tmp.col(j)*xi_l(j);
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
            // update linear alpha star
            // List
            List uasl = update_alphaC(y, sigma, tau_star_l(j,k), gamma_star_l(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_l(j,k));
            alpha_star_l(j,k) = uasl[0];
            int alpha_acc = uasl[1];
            alpha_star_l_acc(j,k) = alpha_star_l_acc(j,k) + alpha_acc;
            m_star_l(j,k) = update_mCsca(xi_star_l(j,k));
            // update xi star linear
            // proposed xi star
            xi_star_bar = xi_star_l(j,k) + mht(1)*randn(); // generate_normal(0, mht(1));
            // compute the new beta_l (temp element)
            omega_l_tmp = omega_l;
            omega_l_tmp(j,k) =  alpha_star_l(j,k) * xi_star_bar;
            alpha_l_tmp.col(j) = alpha_0_l(j)*vecOnes;
            for (int kn = (j+1); kn<p; kn++) {
              alpha_l_tmp.col(j) = alpha_l_tmp.col(j) + X_l.col(kn)*omega_l_tmp(j,kn);
            }
            beta_l_tmp.col(j) = alpha_l_tmp.col(j)*xi_l(j);
            eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
            // update xi linear star
            List uxsl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_star_l(j,k), xi_star_bar, xi_star_l(j,k));
            xi_star_l(j,k) = uxsl[0];
            int xiSacc = uxsl[1];
            xi_star_l_acc(j,k) = xi_star_l_acc(j,k) + xiSacc;
          } else {
            gamma_star_l(j,k) = hyperpar4;
          }
        }
        
        if (ha == 2) {
          if ((gamma_0_l(j) != hyperpar4) && (gamma_0_l(k) != hyperpar4)) {
            // update gamma inclusion parameters
            gamma_star_l(j,k) = update_gammaScaC(pi_star_l, hyperpar4, alpha_star_l(j,k), tau_star_l(j,k));
            // update tau variance
            tau_star_l(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_l(j,k), gamma_star_l(j,k));
            // update alpha star linear and non linear with a MH
            // proposed alpha
            alpha_star_bar = alpha_star_l(j,k) + mht(0)*randn(); // generate_normal(0, mht(0));
            // compute the new beta_l (temp element)
            omega_l_tmp = omega_l;
            omega_l_tmp(j,k) = xi_star_l(j,k) * alpha_star_bar;
            // change the j,k
            alpha_l_tmp.col(j) = alpha_0_l(j)*vecOnes;
            for (int kn = (j+1); kn<p; kn++) {
              alpha_l_tmp.col(j) = alpha_l_tmp.col(j) + X_l.col(kn)*omega_l_tmp(j,kn);
            }
            beta_l_tmp = beta_l;
            // compute linear beta
            beta_l_tmp.col(j) = alpha_l_tmp.col(j)*xi_l(j);
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
            // update linear alpha star
            // List
            List uasl = update_alphaC(y, sigma, tau_star_l(j,k), gamma_star_l(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_l(j,k));
            alpha_star_l(j,k) = uasl[0];
            int alpha_acc = uasl[1];
            alpha_star_l_acc(j,k) = alpha_star_l_acc(j,k) + alpha_acc;
            m_star_l(j,k) = update_mCsca(xi_star_l(j,k));
            // update xi star linear
            // proposed xi star
            xi_star_bar = xi_star_l(j,k) + mht(1)*randn(); // generate_normal(0, mht(1));
            // compute the new beta_l (temp element)
            omega_l_tmp = omega_l;
            omega_l_tmp(j,k) =  alpha_star_l(j,k) * xi_star_bar;
            alpha_l_tmp.col(j) = alpha_0_l(j)*vecOnes;
            for (int kn = (j+1); kn<p; kn++) {
              alpha_l_tmp.col(j) = alpha_l_tmp.col(j) + X_l.col(kn)*omega_l_tmp(j,kn);
            }
            beta_l_tmp.col(j) = alpha_l_tmp.col(j)*xi_l(j);
            eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
            // update xi linear star
            List uxsl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_star_l(j,k), xi_star_bar, xi_star_l(j,k));
            xi_star_l(j,k) = uxsl[0];
            int xiSacc = uxsl[1];
            xi_star_l_acc(j,k) = xi_star_l_acc(j,k) + xiSacc;
          } else {
            gamma_star_l(j,k) = hyperpar4;
          }
        }
        
        if (ha == 0) {
          // update gamma inclusion parameters
          gamma_star_l(j,k) = update_gammaScaC(pi_star_l, hyperpar4, alpha_star_l(j,k), tau_star_l(j,k));
          // update tau variance
          tau_star_l(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_l(j,k), gamma_star_l(j,k));
          // update alpha star linear and non linear with a MH
          // proposed alpha
          alpha_star_bar = alpha_star_l(j,k) + mht(0)*randn(); // generate_normal(0, mht(0));
          // compute the new beta_l (temp element)
          omega_l_tmp = omega_l;
          omega_l_tmp(j,k) = xi_star_l(j,k) * alpha_star_bar;
          // change the j,k
          alpha_l_tmp.col(j) = alpha_0_l(j)*vecOnes;
          for (int kn = (j+1); kn<p; kn++) {
            alpha_l_tmp.col(j) = alpha_l_tmp.col(j) + X_l.col(kn)*omega_l_tmp(j,kn);
          }
          beta_l_tmp = beta_l;
          // compute linear beta
          beta_l_tmp.col(j) = alpha_l_tmp.col(j)*xi_l(j);
          // new linear predictor WITH the proposed alpha star
          eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
          // update linear alpha star
          // List
          List uasl = update_alphaC(y, sigma, tau_star_l(j,k), gamma_star_l(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_l(j,k));
          alpha_star_l(j,k) = uasl[0];
          int alpha_acc = uasl[1];
          alpha_star_l_acc(j,k) = alpha_star_l_acc(j,k) + alpha_acc;
          m_star_l(j,k) = update_mCsca(xi_star_l(j,k));
          // update xi star linear
          // proposed xi star
          xi_star_bar = xi_star_l(j,k) + mht(1)*randn(); // generate_normal(0, mht(1));
          // compute the new beta_l (temp element)
          omega_l_tmp = omega_l;
          omega_l_tmp(j,k) =  alpha_star_l(j,k) * xi_star_bar;
          alpha_l_tmp.col(j) = alpha_0_l(j)*vecOnes;
          for (int kn = (j+1); kn<p; kn++) {
            alpha_l_tmp.col(j) = alpha_l_tmp.col(j) + X_l.col(kn)*omega_l_tmp(j,kn);
          }
          beta_l_tmp.col(j) = alpha_l_tmp.col(j)*xi_l(j);
          eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
          // update xi linear star
          List uxsl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_star_l(j,k), xi_star_bar, xi_star_l(j,k));
          xi_star_l(j,k) = uxsl[0];
          int xiSacc = uxsl[1];
          xi_star_l_acc(j,k) = xi_star_l_acc(j,k) + xiSacc;
        }
      } // end linear k
    } // end linear j
    
    // update gamma non linear
    gamma_0_nl = update_gammaVecC(pi_0_nl, hyperpar4, alpha_0_nl, tau_0_nl);
    if (ha == 2) {
      active_j = find(gamma_0_nl == hyperpar4);
      for (uword j : active_j) {
        if (j < nlp - 1) {
          gamma_star_nl.submat(j, j+1, j, nlp-1).fill(hyperpar4);
        }
      }
    } else if (ha == 1) {
      active_pairs = conv_to<umat>::from(
        (gamma_0_nl * gamma_0_nl.t()) == hyperpar4_sq
      );
      active_pairs = trimatu(active_pairs, 1);
      gamma_star_nl.elem(find(active_pairs > 0)).fill(hyperpar4);
    }
    
    //////////////////// effect modifiers non-linear peNMIG ////////////////////
    for (int j = 0; j<nlp; j++) {
      for (int k = (j+1); k<nlp; k++) {
        // update all the non linear interaction
        if (ha == 1) {
          if ((gamma_0_l(j) != hyperpar4) || (gamma_0_l(k) != hyperpar4)) {
            // update gamma star non linear
            gamma_star_nl(j,k) = update_gammaScaC(pi_star_nl, hyperpar4, alpha_star_nl(j,k), tau_star_nl(j,k));
            // update tau star non linear
            tau_star_nl(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_nl(j,k), gamma_star_nl(j,k));
            // update alpha non linear star
            // proposed alpha
            alpha_star_bar = alpha_star_nl(j,k) + mht(2)*randn(); // generate_normal(0, mht(2));
            omega_nl_tmp = omega_nl;
            omega_nl_tmp(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1)) * alpha_star_bar;
            alpha_nl_tmp.col(j) = alpha_0_nl(j)*vecOnes;
            for (int kn = (j+1); kn<nlp; kn++) {
              alpha_nl_tmp.col(j) = alpha_nl_tmp.col(j) + X_nl.cols(span(cd[kn], cd[kn+1]-1))*omega_nl_tmp(j, span(cd[kn], cd[kn+1]-1)).t();
            }
            // compute the beta non linear temp
            beta_nl_tmp = beta_nl;
            beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl_tmp.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
            // update alpha star non linear
            List uasnl = update_alphaC(y, sigma, tau_star_nl(j,k), gamma_star_nl(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_nl(j,k));
            alpha_star_nl(j,k) = uasnl[0];
            int acc_anl = uasnl[1];
            alpha_star_nl_acc(j,k) = alpha_star_nl_acc(j,k) + acc_anl;
            // update m star non linear vector
            m_star_nl(j, span(cd[k], cd[k+1]-1)) = update_mCvec(xi_star_nl(j, span(cd[k], cd[k+1]-1))).t();
            // update xi star non linear
            // proposed xi star non linear
            xi_star_bar_nl = xi_star_nl(j, span(cd[k], cd[k+1]-1)) + mht(3)*randn(); // generate_normal(0, mht(3));
            // compute the new beta 
            omega_nl_tmp = omega_nl;
            omega_nl_tmp(j, span(cd[k], cd[k+1]-1)) = alpha_star_nl(j,k)*xi_star_bar_nl;
            alpha_nl_tmp.col(j) = alpha_0_nl(j)*vecOnes;
            for (int kn = (j+1); kn<nlp; kn++) {
              alpha_nl_tmp.col(j) = alpha_nl_tmp.col(j) + X_nl.cols(span(cd[kn], cd[kn+1]-1))*omega_nl_tmp(j, span(cd[kn], cd[kn+1]-1)).t();
            }
            // compute the beta non linear temp
            beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl_tmp.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
            // compute the linear predictor this the proposed xi
            eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
            // update xi star non-linear
            List uxsnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_star_nl(j, span(cd[k], cd[k+1]-1)).t(), xi_star_bar_nl.t(), xi_star_nl(j, span(cd[k], cd[k+1]-1)).t());
            rowvec resultXiSnl = uxsnl[0];
            xi_star_nl(j, span(cd[k], cd[k+1]-1)) = resultXiSnl;
            rowvec accxisnl = uxsnl[1];
            xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) = xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) + accxisnl;
          } else {
            gamma_star_nl(j,k) = hyperpar4;
          }
        }
        
        if (ha == 2) {
          if ((gamma_0_l(j) != hyperpar4) && (gamma_0_l(k) != hyperpar4)) {
            // update gamma star non linear
            gamma_star_nl(j,k) = update_gammaScaC(pi_star_nl, hyperpar4, alpha_star_nl(j,k), tau_star_nl(j,k));
            // update tau star non linear
            tau_star_nl(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_nl(j,k), gamma_star_nl(j,k));
            // update alpha non linear star
            // proposed alpha
            alpha_star_bar = alpha_star_nl(j,k) + mht(2)*randn(); // generate_normal(0, mht(2));
            omega_nl_tmp = omega_nl;
            omega_nl_tmp(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1)) * alpha_star_bar;
            alpha_nl_tmp.col(j) = alpha_0_nl(j)*vecOnes;
            for (int kn = (j+1); kn<nlp; kn++) {
              alpha_nl_tmp.col(j) = alpha_nl_tmp.col(j) + X_nl.cols(span(cd[kn], cd[kn+1]-1))*omega_nl_tmp(j, span(cd[kn], cd[kn+1]-1)).t();
            }
            // compute the beta non linear temp
            beta_nl_tmp = beta_nl;
            beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl_tmp.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
            // new linear predictor WITH the proposed alpha star
            eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
            // update alpha star non linear
            List uasnl = update_alphaC(y, sigma, tau_star_nl(j,k), gamma_star_nl(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_nl(j,k));
            alpha_star_nl(j,k) = uasnl[0];
            int acc_anl = uasnl[1];
            alpha_star_nl_acc(j,k) = alpha_star_nl_acc(j,k) + acc_anl;
            // update m star non linear vector
            m_star_nl(j, span(cd[k], cd[k+1]-1)) = update_mCvec(xi_star_nl(j, span(cd[k], cd[k+1]-1))).t();
            // update xi star non linear
            // proposed xi star non linear
            xi_star_bar_nl = xi_star_nl(j, span(cd[k], cd[k+1]-1)) + mht(3)*randn(); // generate_normal(0, mht(3));
            // compute the new beta 
            omega_nl_tmp = omega_nl;
            omega_nl_tmp(j, span(cd[k], cd[k+1]-1)) = alpha_star_nl(j,k)*xi_star_bar_nl;
            alpha_nl_tmp.col(j) = alpha_0_nl(j)*vecOnes;
            for (int kn = (j+1); kn<nlp; kn++) {
              alpha_nl_tmp.col(j) = alpha_nl_tmp.col(j) + X_nl.cols(span(cd[kn], cd[kn+1]-1))*omega_nl_tmp(j, span(cd[kn], cd[kn+1]-1)).t();
            }
            // compute the beta non linear temp
            beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl_tmp.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
            // compute the linear predictor this the proposed xi
            eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
            // update xi star non-linear
            List uxsnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_star_nl(j, span(cd[k], cd[k+1]-1)).t(), xi_star_bar_nl.t(), xi_star_nl(j, span(cd[k], cd[k+1]-1)).t());
            rowvec resultXiSnl = uxsnl[0];
            xi_star_nl(j, span(cd[k], cd[k+1]-1)) = resultXiSnl;
            rowvec accxisnl = uxsnl[1];
            xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) = xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) + accxisnl;
          } else {
            gamma_star_nl(j,k) = hyperpar4;
          }
        }
        
        if (ha == 0) {
          // update gamma star non linear
          gamma_star_nl(j,k) = update_gammaScaC(pi_star_nl, hyperpar4, alpha_star_nl(j,k), tau_star_nl(j,k));
          // update tau star non linear
          tau_star_nl(j,k) = update_tauC(hyperpar(0), hyperpar(1), alpha_star_nl(j,k), gamma_star_nl(j,k));
          // update alpha non linear star
          // proposed alpha
          alpha_star_bar = alpha_star_nl(j,k) + mht(2)*randn(); // generate_normal(0, mht(2));
          omega_nl_tmp = omega_nl;
          omega_nl_tmp(j, span(cd[k], cd[k+1]-1)) = xi_star_nl(j, span(cd[k], cd[k+1]-1)) * alpha_star_bar;
          alpha_nl_tmp.col(j) = alpha_0_nl(j)*vecOnes;
          for (int kn = (j+1); kn<nlp; kn++) {
            alpha_nl_tmp.col(j) = alpha_nl_tmp.col(j) + X_nl.cols(span(cd[kn], cd[kn+1]-1))*omega_nl_tmp(j, span(cd[kn], cd[kn+1]-1)).t();
          }
          // compute the beta non linear temp
          beta_nl_tmp = beta_nl;
          beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl_tmp.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
          // new linear predictor WITH the proposed alpha star
          eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
          // update alpha star non linear
          List uasnl = update_alphaC(y, sigma, tau_star_nl(j,k), gamma_star_nl(j,k), eta_pl_tmp, eta_pl, alpha_star_bar, alpha_star_nl(j,k));
          alpha_star_nl(j,k) = uasnl[0];
          int acc_anl = uasnl[1];
          alpha_star_nl_acc(j,k) = alpha_star_nl_acc(j,k) + acc_anl;
          // update m star non linear vector
          m_star_nl(j, span(cd[k], cd[k+1]-1)) = update_mCvec(xi_star_nl(j, span(cd[k], cd[k+1]-1))).t();
          // update xi star non linear
          // proposed xi star non linear
          xi_star_bar_nl = xi_star_nl(j, span(cd[k], cd[k+1]-1)) + mht(3)*randn(); // generate_normal(0, mht(3));
          // compute the new beta 
          omega_nl_tmp = omega_nl;
          omega_nl_tmp(j, span(cd[k], cd[k+1]-1)) = alpha_star_nl(j,k)*xi_star_bar_nl;
          alpha_nl_tmp.col(j) = alpha_0_nl(j)*vecOnes;
          for (int kn = (j+1); kn<nlp; kn++) {
            alpha_nl_tmp.col(j) = alpha_nl_tmp.col(j) + X_nl.cols(span(cd[kn], cd[kn+1]-1))*omega_nl_tmp(j, span(cd[kn], cd[kn+1]-1)).t();
          }
          // compute the beta non linear temp
          beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl_tmp.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
          // compute the linear predictor this the proposed xi
          eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
          // update xi star non-linear
          List uxsnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_star_nl(j, span(cd[k], cd[k+1]-1)).t(), xi_star_bar_nl.t(), xi_star_nl(j, span(cd[k], cd[k+1]-1)).t());
          rowvec resultXiSnl = uxsnl[0];
          xi_star_nl(j, span(cd[k], cd[k+1]-1)) = resultXiSnl;
          rowvec accxisnl = uxsnl[1];
          xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) = xi_star_nl_acc(j, span(cd[k], cd[k+1]-1)) + accxisnl;
        }
      } // end non linear k
    } // end non linear j
    
    
    // rescaling linear
    for (int j = 0; j<p; j++) {
      for (int k = (j+1); k<p; k++) {
        // linear
        sFct = 1/accu(abs(xi_star_l(j,k)));
        xi_star_l(j,k) = xi_star_l(j,k)*sFct;
        alpha_star_l(j,k) = alpha_star_l(j,k)/sFct;
        // compute linear omega
        omega_l(j,k) = alpha_star_l(j,k)*xi_star_l(j,k);
      }
    }
    
    // rescaling non linear
    for (int j = 0; j<nlp; j++) {
      for (int k = (j+1); k<nlp; k++) {
        // non linear
        sFct = d(k)/accu(abs(xi_star_nl(j, span(cd[k], cd[k+1]-1))));
        xi_star_nl(j, span(cd[k], cd[k+1]-1)) = sFct*xi_star_nl(j, span(cd[k], cd[k+1]-1));
        alpha_star_nl(j,k) = alpha_star_nl(j,k)/sFct;
        // compute non linear omega
        omega_nl(j, span(cd[k], cd[k+1]-1)) = alpha_star_nl(j,k)*xi_star_nl(j, span(cd[k], cd[k+1]-1));
      }
    }
    
    // update pi 0 linear
    pi_0_l = update_piNSC(hyperpar(5), hyperpar(6), gamma_0_l, hyperpar4);
    // update pi 0 non linear
    pi_0_nl = update_piNSC(hyperpar(7), hyperpar(8), gamma_0_nl, hyperpar4);
    // update gamma 0 linear
    
    // If ha = 2 and I remove the main effect j, I must also eliminate all interactions 
    // involving it, including gamma.
    // If ha = 1 and I remove the main effect j, I must ensure that the other main effect 
    // is present; otherwise, I remove the interaction.
    // gamma_0_l = update_gammaVecC(pi_0_l, hyperpar4, alpha_0_l, tau_0_l);
    // update gamma 0 non linear
    // gamma_0_nl = update_gammaVecC(pi_0_nl, hyperpar4, alpha_0_nl, tau_0_nl);
    
    // update tau 0 linear
    for (int j = 0; j<p; j++) {
      tau_0_l(j) = update_tauC(hyperpar(0), hyperpar(1), alpha_0_l(j), gamma_0_l(j));
    }
    // update tau 0 non linear
    for (int j = 0; j<nlp; j++) {
      tau_0_nl(j) = update_tauC(hyperpar(0), hyperpar(1), alpha_0_nl(j), gamma_0_nl(j));
    }
    
    // update alpha 0 linear
    for (int j = 0; j<p; j++) {
      // proposed alpha
      alpha_0_bar = alpha_0_l(j) + mht(4)*randn(); // generate_normal(0, mht(4));
      alpha_l_tmp.col(j) = alpha_0_bar*vecOnes;
      for (int k = (j+1); k<p; k++) {
        alpha_l_tmp.col(j) = alpha_l_tmp.col(j) + X_l.col(k)*omega_l(j, k);
      }
      // compute linear beta
      beta_l_tmp = beta_l;
      beta_l_tmp.col(j) = alpha_l_tmp.col(j)*xi_l(j);
      eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
      // update alpha 0 linear
      List ua0l = update_alphaC(y, sigma, tau_0_l(j), gamma_0_l(j), eta_pl_tmp, eta_pl, alpha_0_bar, alpha_0_l(j));
      alpha_0_l(j) = ua0l[0];
      int accAlpha0 = ua0l[1];
      alpha_0_l_acc(j) = alpha_0_l_acc(j) + accAlpha0;
    }
    
    // update alpha 0 non linear
    for (int j = 0; j<nlp; j++) {
      // proposed alpha
      alpha_0_bar = alpha_0_nl(j) + mht(5)*randn(); // generate_normal(0, mht(5));
      // alpha_nl_tmp = alpha_nl;
      alpha_nl_tmp.col(j) = alpha_0_bar*vecOnes;
      for (int k = (j+1); k<nlp; k++) {
        alpha_nl_tmp.col(j) = alpha_nl_tmp.col(j) + X_nl.cols(span(cd[k], cd[k+1]-1))*omega_nl(j, span(cd[k], cd[k+1]-1)).t();
      }
      // compute non linear beta
      beta_nl_tmp = beta_nl;
      beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl_tmp.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
      eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
      List ua0nl = update_alphaC(y, sigma, tau_0_nl(j), gamma_0_nl(j), eta_pl_tmp, eta_pl, alpha_0_bar, alpha_0_nl(j));
      alpha_0_nl(j) = ua0nl[0];
      int accAnlpha0 = ua0nl[1];
      alpha_0_nl_acc(j) = alpha_0_nl_acc(j) + accAnlpha0;
    }
    
    // compute alpha linear
    for (int j = 0; j<p; j++) {
      alpha_l.col(j) = alpha_0_l(j)*vecOnes;
      for (int k = (j+1); k<p; k++) {
        alpha_l.col(j) = alpha_l.col(j) + X_l.col(k)*omega_l(j,k);
      }
    }
    
    // compute alpha non linear
    for (int j = 0; j<nlp; j++) {
      alpha_nl.col(j) = alpha_0_nl(j)*vecOnes;
      for (int k = (j+1); k<nlp; k++) {
        alpha_nl.col(j) = alpha_nl.col(j) + X_nl.cols(span(cd[k], cd[k+1]-1))*omega_nl(j, span(cd[k], cd[k+1]-1)).t();
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
      xi_star = xi_l(j) + mht(6)*randn(); // generate_normal(0, mht(6));
      // compute the new beta linear
      beta_l_tmp = beta_l; 
      beta_l_tmp.col(j) = alpha_l.col(j)*xi_star;
      // compute the new linear predictor with the proposed xi
      eta_pl_tmp = eta_pl + X_l.col(j)%(beta_l_tmp.col(j) - beta_l.col(j));
      // update xi linear
      List uxl = update_xiLC(y, eta_pl_tmp, eta_pl, sigma, m_l(j), xi_star, xi_l(j));
      xi_l(j) = uxl[0];
      int acc_xil = uxl[1];
      xi_l_acc(j) = xi_l_acc(j) + acc_xil;
    }
    
    // update xi non linear
    xi_starnl = xi_nl + as<vec>(wrap(Rcpp::rnorm(q, 0, mht(7))));
    for (int j = 0; j<nlp; j++) {
      beta_nl_tmp = beta_nl; 
      beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) = alpha_nl.col(j)*xi_starnl(span(cd[j], cd[j+1]-1)).t();
      eta_pl_tmp = eta_pl + sum(X_nl.cols(span(cd[j], cd[j+1]-1))%(beta_nl_tmp.cols(span(cd[j], cd[j+1]-1)) - beta_nl.cols(span(cd[j], cd[j+1]-1))), 1);
      List uxnl = update_xiNLC(y, eta_pl_tmp, eta_pl, sigma, m_nl(span(cd[j], cd[j+1]-1)), xi_starnl(span(cd[j], cd[j+1]-1)), xi_nl(span(cd[j], cd[j+1]-1)));
      vec resXnl = uxnl[0];
      xi_nl(span(cd[j], cd[j+1]-1)) = resXnl;
      vec accXinl = uxnl[1];
      xi_nl_acc(span(cd[j], cd[j+1]-1)) = xi_nl_acc(span(cd[j], cd[j+1]-1)) + accXinl;
    }
    
    // rescale alpha and xi linear
    for (int j = 0; j<p; j++) {
      // scaling factor
      sFct = 1/accu(abs(xi_l(j)));
      // rescale linear xi
      xi_l(j) = xi_l(j)*sFct;
      // rescale alpha linear
      alpha_l.col(j) = alpha_l.col(j)/sFct;
    }
    
    // rescale alpha and xi non linear
    for (int j = 0; j<nlp; j++) {
      // scaling factor
      sFct = d(j)/accu(abs(xi_nl(span(cd[j], cd[j+1]-1))));
      // rescale non linear xi
      xi_nl(span(cd[j], cd[j+1]-1)) = sFct*xi_nl(span(cd[j], cd[j+1]-1));
      // rescale non linear alpha
      alpha_nl.col(j) = alpha_nl.col(j)/sFct;
    }
    
    // update beta linear and non linear after rescaling alpha and xi
    // beta linear and non linear
    for (int j = 0; j<p; j++) {
      beta_l.col(j) = alpha_l.col(j)*xi_l(j);
    }
    
    // non linear
    for (int j = 0; j<nlp; j++) {
      beta_nl.cols(span(cd[j], cd[j+1]-1)) = alpha_nl.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
    }
    
    // compute the linear predictor
    eta_pl = compLinPred(nobs, eta0, X_l, beta_l, X_nl, beta_nl);
    // update intercept
    // linear predictor without intercept
    vec eta_noInt = eta_pl - eta0*vecOnes;
    eta0 = updateInterceptC(y, nobs, eta_noInt, sigma);
    // update linear predictor
    eta_pl = eta0*vecOnes + eta_noInt;
    // update sigma variance
    sigma = update_sigmaC(y, eta_pl, hyperpar(2), hyperpar(3), nobs);
    // log-likelihood
    double logLik =  accu(dnormLogVec(y, eta_pl, sqrt(sigma)));
    // store resutls
    if(t%thin == 0 && t > burnin-1) { // we start from 0
      if (detail == true) {
        PI_S_l(idx) = pi_star_l;
        PI_S_nl(idx) = pi_star_nl;
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
        PI_0_l(idx) = pi_0_l;
        PI_0_nl(idx) = pi_0_nl;
        ETA_PL.row(idx) = eta_pl.t();
        ETA0(idx) = eta0;
        SIGMA(idx) = sigma;
        LOGLIKELIHOOD(idx) = logLik;
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
      } else {
        LOGLIKELIHOOD(idx) = logLik;
        SIGMA(idx) = sigma;
        gamma_0_l_m = gamma_0_l_m + gamma_0_l;
        gamma_0_nl_m = gamma_0_nl_m + gamma_0_nl;
        gamma_star_l_m = gamma_star_l_m + gamma_star_l;
        gamma_star_nl_m = gamma_star_nl_m + gamma_star_nl;
        eta_pl_m = eta_pl_m + eta_pl;
      }
      if (n_val != 0) {
        // compute alpha linear
        for (int j = 0; j<p; j++) {
          alpha_val_l.col(j) = alpha_0_l(j)*vecOnesVal;
          for (int k = (j+1); k<p; k++) {
            alpha_val_l.col(j) = alpha_val_l.col(j) + X_val_l.col(k)*omega_l(j,k);
          }
        }
        // compute alpha non linear
        for (int j = 0; j<nlp; j++) {
          alpha_val_nl.col(j) = alpha_0_nl(j)*vecOnesVal;
          for (int k = (j+1); k<nlp; k++) {
            alpha_val_nl.col(j) = alpha_val_nl.col(j) + X_val_nl.cols(span(cd[k], cd[k+1]-1))*omega_nl(j, span(cd[k], cd[k+1]-1)).t();
          }
        }
        // beta linear and non linear
        for (int j = 0; j<p; j++) {
          beta_val_l.col(j) = alpha_val_l.col(j)*xi_l(j);
        }
        // non linear
        for (int j = 0; j<nlp; j++) {
          beta_val_nl.cols(span(cd[j], cd[j+1]-1)) = alpha_val_nl.col(j)*xi_nl(span(cd[j], cd[j+1]-1)).t();
        }
        // compute linear predictor
        eta_pl_val = compLinPred(n_val, eta0, X_val_l, beta_val_l, X_val_nl, beta_val_nl);
        // y_tilde
        for (int i = 0; i<n_val; i++) {
          y_tilde(i) = generate_normal(eta_pl_val(i), sqrt(sigma));
        }
        Y_TILDE.row(idx) = y_tilde.t();
        y_tilde_m = y_tilde_m + y_tilde;
      }
      idx = idx + 1;
    }
    if (pb == true && t%500 == 0 && t > 0) {
      std::cout << "Iteration: " << t << " (of " << iter << ")\n";
    }
  } // end iteration
  if (pb == true) {
    std::cout << "End MCMC!\n";
  }
  ////////////////////////////////////////////////////
  //////////////////// End MCMC //////////////////////
  ///////////////////////////////////////////////////
  // Time 
  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  // mean no detail
  gamma_0_l_m = gamma_0_l_m/nout;
  gamma_0_nl_m = gamma_0_nl_m/nout;
  gamma_star_l_m = gamma_star_l_m/nout;
  gamma_star_nl_m = gamma_star_nl_m/nout;
  eta_pl_m = eta_pl_m/nout;
  y_tilde_m = y_tilde_m/nout;
  // std::cout << "The computational time for the entire MCMC is " << duration/1000000;
  Rcpp::List acc_list = List::create(Named("acc_a_s_l") = alpha_star_l_acc/iter,
                                     Named("acc_a_s_nl") = alpha_star_nl_acc/iter,
                                     Named("acc_xi_s_l") = xi_star_l_acc/iter,
                                     Named("acc_xi_s_nl") = xi_star_nl_acc/iter,
                                     Named("acc_a_0_l") = alpha_0_l_acc/iter, 
                                     Named("acc_a_0_nl") = alpha_0_nl_acc/iter, 
                                     Named("acc_xi_l") = xi_l_acc/iter, 
                                     Named("acc_xi_nl") = xi_nl_acc/iter);
  if (detail == true) {
    return List::create(Named("d") = d,
                        Named("X_lin") = X_l,
                        Named("X_nl") = X_nl,
                        // Named("intercept") = ETA0,
                        Named("linear_predictor") = ETA_PL,
                        // Named("Beta_l") = BETA_l,
                        // Named("Beta_nl") = BETA_nl,
                        Named("alpha_0_l") = ALPHA_0_l,
                        Named("alpha_0_nl") = ALPHA_0_nl,
                        // Named("alpha_star_l") = ALPHA_S_l,
                        // Named("alpha_star_nl") = ALPHA_S_nl,
                        // Named("m_l") = M_l,
                        // Named("m_nl") = M_nl,
                        // Named("m_star_l") = M_S_l,
                        // Named("m_star_nl") = M_S_nl,
                        Named("xi_l") = XI_l,
                        Named("xi_nl") = XI_nl,
                        Named("omega_l") = OMEGA_l,
                        Named("omega_nl") = OMEGA_nl,
                        // Named("xi_star_l") = XI_S_l,
                        // Named("xi_star_nl") = XI_S_nl,
                        // Named("tau_0_l") = TAU_0_l,
                        // Named("tau_0_nl") = TAU_0_nl,
                        // Named("tau_star_l") = TAU_S_l,
                        // Named("tau_star_nl") = TAU_S_nl
                        Named("gamma_0_l") = GAMMA_0_l,
                        Named("gamma_0_nl") = GAMMA_0_nl,
                        Named("gamma_star_l") = GAMMA_S_l,
                        Named("gamma_star_nl") = GAMMA_S_nl,
                        // Named("pi_0_l") = PI_0_l,
                        // Named("pi_0_nl") = PI_0_nl,
                        // Named("pi_star_l") = PI_S_l,
                        // Named("pi_star_nl") = PI_S_nl,
                        Named("sigma") = SIGMA,
                        Named("LogLikelihood") = LOGLIKELIHOOD,
                        // Named("acc_rate") = acc_list,
                        Named("y_oos") = Y_TILDE,
                        Named("Execution_Time") = duration/1000000
    );
  } else {
    return List::create(Named("linear_predictor") = eta_pl_m,
                        Named("gamma_0_l") = gamma_0_l_m,
                        Named("gamma_0_nl") = gamma_0_nl_m,
                        Named("gamma_star_l") = gamma_star_l_m,
                        Named("gamma_star_nl") = gamma_star_nl_m,
                        Named("LogLikelihood") = LOGLIKELIHOOD,
                        Named("y_oos") = y_tilde_m,
                        Named("y_tilde") = Y_TILDE,
                        Named("sigma") = SIGMA,
                        Named("acc_rate") = acc_list,
                        Named("Execution_Time") = duration/1000000
    );
  }
}

