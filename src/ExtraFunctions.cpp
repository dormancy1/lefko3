// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;


// Index of Functions
// 1. NumericVector ricker3  Two-Parameter Ricker Function
// 2. NumericVector beverton3  Two-Parameter Beverton-Holt Function
// 3. NumericVector usher3  Two-Parameter Usher Function
// 4. NumericVector logistic3  Two-Parameter logistic Function
// 5. int .integer_test  Function to Test Whether a Numeric Vector Is Composed Only of Integers
// 6. int .binomial_test  Function to Test Whether an Integer Vector Is Entirely Binomial


//' Two-Parameter Ricker Function
//' 
//' Function \code{ricker3()} creates a vector of values produced by the two-
//' parameter Ricker function as applied with a user-specified time lag. The
//' two-parameter Ricker function is given as 
//' \eqn{\phi_{t+1} = \phi_t \alpha e^{-\beta n_t}}. Here, if no
//' \code{separate_N} vector is provided, then \eqn{n_t = \phi_t}.
//' 
//' @name ricker3
//' 
//' @param start_value A positive number to start the return vector in time 0.
//' @param alpha The alpha parameter in the two-parameter Ricker function. Must
//' be non-negative.
//' @param beta The beta parameter in the two-parameter Ricker function.
//' @param time_steps The number of time steps to run the projection. Must be a
//' positive integer.
//' @param time_lag A positive integer denoting the number of time steps back
//' for the value of phi in the two-parameter Ricker function.
//' @param pre0_subs A logical value indicating whether to use a number other
//' than that given in \code{start_value} for values of phi lagged from times
//' prior to time 0.
//' @param pre0_value A positive number to use for phi lagged from times prior
//' to time 0. Only used if \code{pre0_subs = TRUE}.
//' @param substoch An integer value indicating the kind of substochasticity to
//' use. Values include: \code{0}, no substochasticity enforced (the default);
//' \code{1}, all numbers must be non-negative; and \code{2}, all numbers should
//' be forced to the interval [0, 1].
//' @param separate_N An optional numeric vector with values of N in each time,
//' if phi is to be treated as different from N in the two-parameter model.
//' 
//' @return A numeric vector of values showing values projected under the two-
//' parameter Ricker function.
//' 
//' @examples
//' trial_run1 <- ricker3(1, alpha = 0.5, beta = -0.009)
//' plot(trial_run1)
//' 
//' trial_run2 <- ricker3(1, alpha = 0.5, beta = 0.009)
//' plot(trial_run2)
//' 
//' trial_run3 <- ricker3(1, alpha = 1, beta = -0.009)
//' plot(trial_run3)
//' 
//' trial_run4 <- ricker3(1, alpha = 1, beta = 0.009)
//' plot(trial_run4)
//' 
//' trial_run5 <- ricker3(1, alpha = 5, beta = -0.009)
//' plot(trial_run5)
//' 
//' trial_run6 <- ricker3(1, alpha = 5, beta = 0.009)
//' plot(trial_run6)
//' 
//' used_Ns <- c(10, 15, 12, 14, 14, 150, 15, 1, 5, 7, 9, 14, 13, 16, 17, 19,
//'   25, 26)
//' trial_run7 <- ricker3(1, alpha = 1, beta = -0.009, separate_N = used_Ns)
//' plot(trial_run7)
//' 
//' @export ricker3
// [[Rcpp::export(ricker3)]]
Rcpp::NumericVector ricker3(double start_value, double alpha, double beta,
  int time_steps = 100, int time_lag = 1, bool pre0_subs = false,
  double pre0_value = 0.0, int substoch = 0,
  Nullable<NumericVector> separate_N = R_NilValue) {
  
  NumericVector sepN;
  
  bool base_N = false;
  
  if (start_value <= 0.0) {
    throw Rcpp::exception("Option start_value must be positive.", false);
  }
  if (alpha < 0.0) {
    throw Rcpp::exception("Option alpha must be non-negative.", false);
  }
  if (time_lag < 1) {
    throw Rcpp::exception("Option time_lag must be positive.", false);
  }
  if (pre0_subs && pre0_value <= 0.0) {
    throw Rcpp::exception("Option pre0_value must be positive if pre0_subs is set to TRUE",
      false);
  }
  if (substoch < 0 || substoch > 2) {
    throw Rcpp::exception("Option substoch must equal 0, 1, or 2", false);
  }
  
  if (separate_N.isNotNull()) {
    base_N = true;
    NumericVector sepN_temp(separate_N);
    sepN = sepN_temp;
    
    int sepN_length = sepN.length();
    if (sepN_length != (time_steps - 1)) {
      time_steps = sepN_length - 1;
      Rf_warningcall(R_NilValue,
        "Resetting time_steps to length of separate_N - 1.");
    }
  }
  
  if (time_steps < 1) {
    throw Rcpp::exception("Option time_steps must be positive.", false);
  }
  
  NumericVector output(time_steps + 1);
  double used_phi {0.0};
  double nt {0.0};
  
  output(0) = start_value;
  
  for (int i = 1; i <= time_steps; i++ ) {
    used_phi = 0.0;
    nt = 0.0;
    
    if ((i - time_lag) < 0) {
      if (pre0_subs) {
        used_phi = pre0_value;
      } else {
        used_phi = start_value;
      }
      
      if (base_N) {
        nt = sepN(0);
      } else {
        nt = used_phi;
      }
    } else {
      used_phi = output(i - time_lag);
      
      if (base_N) {
        nt = sepN(i - time_lag);
      } else {
        nt = used_phi;
      }
    }
    
    output(i) = used_phi * alpha * exp((-1 * beta * nt));
    
    if (substoch > 0) {
      if (output(i) < 0.0) {
        output(i) = 0.0;
      } else if (output(i) > 1.0 && substoch == 2) {
        output(i) = 1.0;
      }
    }
  }
  return output;
}

//' Two-Parameter Beverton-Holt Function
//' 
//' Function \code{beverton3()} creates a vector of values produced by the two-
//' parameter Beverton-Holt function as applied with a user-specified time lag.
//' The two-parameter Beverton-Holt function is given as 
//' \eqn{\phi_{t+1} = \phi_t \alpha / (1 + \beta n_t)}. Here, if no
//' \code{separate_N} vector is provided, then \eqn{n_t = \phi_t}.
//' 
//' @name beverton3
//' 
//' @param start_value A positive number to start the return vector in time 0.
//' @param alpha The alpha parameter in the two-parameter Beverton-Holt
//' function. Must be non-negative.
//' @param beta The beta parameter in the two-parameter Beverton-Holt function.
//' Must be non-negative.
//' @param time_steps The number of time steps to run the projection. Must be a
//' positive integer.
//' @param time_lag A positive integer denoting the number of time steps back
//' for the value of phi in the two-parameter Beverton-Holt function.
//' @param pre0_subs A logical value indicating whether to use a number other
//' than that given in \code{start_value} for values of phi lagged from times
//' prior to time 0.
//' @param pre0_value A positive number to use for phi lagged from times prior
//' to time 0. Only used if \code{pre0_subs = TRUE}.
//' @param substoch An integer value indicating the kind of substochasticity to
//' use. Values include: \code{0}, no substochasticity enforced (the default);
//' \code{1}, all numbers must be non-negative; and \code{2}, all numbers should
//' be forced to the interval [0, 1].
//' @param separate_N An optional numeric vector with values of N in each time,
//' if phi is to be treated as different from N in the two-parameter model.
//' 
//' @return A numeric vector of values showing values projected under the two-
//' parameter Beverton-Holt function.
//' 
//' @examples
//' trial_run1 <- beverton3(1, alpha = 0.5, beta = 0.009)
//' plot(trial_run1)
//' 
//' trial_run2 <- beverton3(1, alpha = 0.5, beta = 0.9)
//' plot(trial_run2)
//' 
//' trial_run3 <- beverton3(1, alpha = 1, beta = 0.009)
//' plot(trial_run3)
//' 
//' trial_run4 <- beverton3(1, alpha = 1, beta = 0.9)
//' plot(trial_run4)
//' 
//' trial_run5 <- beverton3(1, alpha = 5, beta = 0.009)
//' plot(trial_run5)
//' 
//' trial_run6 <- beverton3(1, alpha = 5, beta = 0.9)
//' plot(trial_run6)
//' 
//' used_Ns <- c(10, 15, 12, 14, 14, 150, 15, 1, 5, 7, 9, 14, 13, 16, 17, 19,
//'   25, 26)
//' trial_run7 <- beverton3(1, alpha = 1, beta = 0.009, separate_N = used_Ns)
//' plot(trial_run7)
//' 
//' @export beverton3
// [[Rcpp::export(beverton3)]]
Rcpp::NumericVector beverton3(double start_value, double alpha, double beta,
  int time_steps = 100, int time_lag = 1, bool pre0_subs = false,
  double pre0_value = 0.0, int substoch = 0,
  Nullable<NumericVector> separate_N = R_NilValue) {
  
  NumericVector sepN;
  
  bool base_N = false;
  
  if (start_value <= 0.0) {
    throw Rcpp::exception("Option start_value must be positive.", false);
  }
  if (alpha < 0.0) {
    throw Rcpp::exception("Option alpha must be non-negative.", false);
  }
  if (beta < 0.0) {
    throw Rcpp::exception("Option beta must be non-negative.", false);
  }
  if (time_lag < 1) {
    throw Rcpp::exception("Option time_lag must be positive.", false);
  }
  if (pre0_subs && pre0_value <= 0.0) {
    throw Rcpp::exception("Option pre0_value must be positive if pre0_subs is set to TRUE",
      false);
  }
  if (substoch < 0 || substoch > 2) {
    throw Rcpp::exception("Option substoch must equal 0, 1, or 2", false);
  }
  
  if (separate_N.isNotNull()) {
    base_N = true;
    NumericVector sepN_temp(separate_N);
    sepN = sepN_temp;
    
    int sepN_length = sepN.length();
    if (sepN_length != (time_steps - 1)) {
      time_steps = sepN_length - 1;
      Rf_warningcall(R_NilValue,
        "Resetting time_steps to length of separate_N - 1.");
    }
  }
  
  if (time_steps < 1) throw Rcpp::exception("Option time_steps must be positive.",
    false);
  
  NumericVector output(time_steps + 1);
  double used_phi {0.0};
  double nt {0.0};
  
  output(0) = start_value;
  
  for (int i = 1; i <= time_steps; i++ ) {
    used_phi = 0.0;
    nt = 0.0;
    
    if ((i - time_lag) < 0) {
      if (pre0_subs) {
        used_phi = pre0_value;
      } else {
        used_phi = start_value;
      }
      
      if (base_N) {
        nt = sepN(0);
      } else {
        nt = used_phi;
      }
    } else {
      used_phi = output(i - time_lag);
      
      if (base_N) {
        nt = sepN(i - time_lag);
      } else {
        nt = used_phi;
      }
    }
    
    output(i) = used_phi * alpha / (1 + beta * nt);
    
    if (substoch > 0) {
      if (output(i) < 0.0) {
        output(i) = 0.0;
      } else if (output(i) > 1.0 && substoch == 2) {
        output(i) = 1.0;
      }
    }
 }
  return output;
}

//' Two-Parameter Usher Function
//' 
//' Function \code{usher3()} creates a vector of values produced by the two-
//' parameter Usher function as applied with a user-specified time lag.
//' The Usher function is given as 
//' \eqn{\phi_{t+1} = \phi_t / (1 + e^{\alpha n_t + \beta})}. Here,
//' if no \code{separate_N} vector is provided, then \eqn{n_t = \phi_t}.
//' 
//' @name usher3
//' 
//' @param start_value A positive number to start the return vector in time 0.
//' @param alpha The alpha parameter in the two-parameter Usher
//' function.
//' @param beta The beta parameter in the two-parameter Usher function.
//' @param time_steps The number of time steps to run the projection. Must be a
//' positive integer.
//' @param time_lag A positive integer denoting the number of time steps back
//' for the value of phi in the two-parameter Usher function.
//' @param pre0_subs A logical value indicating whether to use a number other
//' than that given in \code{start_value} for values of phi lagged from times
//' prior to time 0.
//' @param pre0_value A positive number to use for phi lagged from times prior
//' to time 0. Only used if \code{pre0_subs = TRUE}.
//' @param substoch An integer value indicating the kind of substochasticity to
//' use. Values include: \code{0}, no substochasticity enforced (the default);
//' \code{1}, all numbers must be non-negative; and \code{2}, all numbers should
//' be forced to the interval [0, 1].
//' @param separate_N An optional numeric vector with values of N in each time,
//' if phi is to be treated as different from N in the two-parameter model.
//' 
//' @return A numeric vector of values showing values projected under the two-
//' parameter Usher function.
//' 
//' @examples
//' trial_run1 <- usher3(1, alpha = -0.5, beta = 0.005)
//' plot(trial_run1)
//' 
//' trial_run2 <- usher3(1, alpha = 0.5, beta = 0.005)
//' plot(trial_run2)
//' 
//' trial_run3 <- usher3(1, alpha = -5, beta = 0.005)
//' plot(trial_run3)
//' 
//' trial_run4 <- usher3(1, alpha = 5, beta = 0.005)
//' plot(trial_run4)
//' 
//' trial_run5 <- usher3(1, alpha = -25, beta = 0.005)
//' plot(trial_run5)
//' 
//' trial_run6 <- usher3(1, alpha = 25, beta = 0.005)
//' plot(trial_run6)
//' 
//' used_Ns <- c(10, 15, 12, 14, 14, 150, 15, 1, 5, 7, 9, 14, 13, 16, 17, 19,
//'   25, 26)
//' trial_run7 <- usher3(1, alpha = -0.5, beta = 0.005, separate_N = used_Ns)
//' plot(trial_run7)
//' 
//' @export usher3
// [[Rcpp::export(usher3)]]
Rcpp::NumericVector usher3(double start_value, double alpha, double beta,
  int time_steps = 100, int time_lag = 1, bool pre0_subs = false,
  double pre0_value = 0.0, int substoch = 0,
  Nullable<NumericVector> separate_N = R_NilValue) {
  
  NumericVector sepN;
  
  bool base_N = false;
  
  if (start_value <= 0.0) {
    throw Rcpp::exception("Option start_value must be positive.", false);
  }
  if (time_lag < 1) {
    throw Rcpp::exception("Option time_lag must be positive.", false);
  }
  if (pre0_subs && pre0_value <= 0.0) {
    throw Rcpp::exception("Option pre0_value must be positive if pre0_subs is set to TRUE",
      false);
  }
  if (substoch < 0 || substoch > 2) {
    throw Rcpp::exception("Option substoch must equal 0, 1, or 2", false);
  }
  
  if (separate_N.isNotNull()) {
    base_N = true;
    NumericVector sepN_temp(separate_N);
    sepN = sepN_temp;
    
    int sepN_length = sepN.length();
    if (sepN_length != (time_steps - 1)) {
      time_steps = sepN_length - 1;
      Rf_warningcall(R_NilValue,
        "Resetting time_steps to length of separate_N - 1.");
    }
  }
  
  if (time_steps < 1) throw Rcpp::exception("Option time_steps must be positive.",
    false);
  
  NumericVector output(time_steps + 1);
  double used_phi {0.0};
  double nt {0.0};
  
  output(0) = start_value;
  
  for (int i = 1; i <= time_steps; i++ ) {
    used_phi = 0.0;
    nt = 0.0;
    
    if ((i - time_lag) < 0) {
      if (pre0_subs) {
        used_phi = pre0_value;
      } else {
        used_phi = start_value;
      }
      
      if (base_N) {
        nt = sepN(0);
      } else {
        nt = used_phi;
      }
    } else {
      used_phi = output(i - time_lag);
      
      if (base_N) {
        nt = sepN(i - time_lag);
      } else {
        nt = used_phi;
      }
    }
    
    output(i) = used_phi / (1 + exp(alpha * nt + beta));
    
    if (substoch > 0) {
      if (output(i) < 0.0) {
        output(i) = 0.0;
      } else if (output(i) > 1.0 && substoch == 2) {
        output(i) = 1.0;
      }
    }
  }
  return output;
}

//' Two-Parameter logistic Function
//' 
//' Function \code{logistic3()} creates a vector of values produced by the
//' logistic function as applied with a user-specified time lag. The logistic
//' function is given as 
//' \eqn{\phi_{t+1} = \phi_t \lambda (1 - n_t / K)}. Here,
//' if no \code{separate_N} vector is provided, then \eqn{n_t = \phi_t}. If
//' \eqn{\lambda} is not provided, then it defaults to \code{1.0}.
//' 
//' @name logistic3
//' 
//' @param start_value A positive number to start the return vector in time 0.
//' @param alpha The carrying capacity K.
//' @param beta If set to some positive number, then this number is the maximum
//' value of phi to enforce. Otherwise, equals \code{0} and enforces no limit.
//' @param lambda The value of the discrete population growth rate to use.
//' Equal to the natural logarithm of the instantaneous growth rate, r.
//' @param time_steps The number of time steps to run the projection. Must be a
//' positive integer.
//' @param time_lag A positive integer denoting the number of time steps back
//' for the value of phi in the logistic function.
//' @param pre0_subs A logical value indicating whether to use a number other
//' than that given in \code{start_value} for values of phi lagged from times
//' prior to time 0.
//' @param pre0_value A positive number to use for phi lagged from times prior
//' to time 0. Only used if \code{pre0_subs = TRUE}.
//' @param substoch An integer value indicating the kind of substochasticity to
//' use. Values include: \code{0}, no substochasticity enforced (the default);
//' \code{1}, all numbers must be non-negative; and \code{2}, all numbers should
//' be forced to the interval [0, 1].
//' @param separate_N An optional numeric vector with values of N in each time,
//' if phi is to be treated as different from N in the logistic model.
//' 
//' @return A numeric vector of values showing values projected under the-
//' logistic function.
//' 
//' @examples
//' trial_run1 <- logistic3(1, alpha = 5)
//' plot(trial_run1)
//' 
//' trial_run2 <- logistic3(1, alpha = 5, beta = 5)
//' plot(trial_run2)
//' 
//' trial_run3 <- logistic3(1, alpha = 100)
//' plot(trial_run3)
//' 
//' trial_run4 <- logistic3(1, alpha = 100, beta = 50)
//' plot(trial_run4)
//' 
//' trial_run5 <- logistic3(1, alpha = 500)
//' plot(trial_run5)
//' 
//' trial_run6 <- logistic3(1, alpha = 500, beta = 501)
//' plot(trial_run6)
//' 
//' used_Ns <- c(10, 15, 12, 14, 14, 150, 15, 1, 5, 7, 9, 14, 13, 16, 17, 19,
//'   25, 26)
//' trial_run7 <- logistic3(1, alpha = 500, beta = 501, separate_N = used_Ns)
//' plot(trial_run7)
//' 
//' @export logistic3
// [[Rcpp::export(logistic3)]]
Rcpp::NumericVector logistic3(double start_value, double alpha,
  double beta = 0.0, double lambda = 1.0, int time_steps = 100, int time_lag = 1,
  bool pre0_subs = false, double pre0_value = 0.0, int substoch = 0,
  Nullable<NumericVector> separate_N = R_NilValue) {
  
  NumericVector sepN;
  
  bool base_N = false;
  
  if (start_value <= 0.0) {
    throw Rcpp::exception("Option start_value must be positive.", false);
  }
  if (alpha <= 0.0) {
    throw Rcpp::exception("Option alpha must be positive.", false);
  }
  if (lambda < 0.0) {
    throw Rcpp::exception("Option lambda must be non-negative.", false);
  }
  if (time_lag < 1) {
    throw Rcpp::exception("Option time_lag must be positive.", false);
  }
  if (pre0_subs && pre0_value <= 0.0) {
    throw Rcpp::exception("Option pre0_value must be positive if pre0_subs is set to TRUE",
      false);
  }
  if (substoch < 0 || substoch > 2) {
    throw Rcpp::exception("Option substoch must equal 0, 1, or 2", false);
  }
  
  if (separate_N.isNotNull()) {
    base_N = true;
    NumericVector sepN_temp(separate_N);
    sepN = sepN_temp;
    
    int sepN_length = sepN.length();
    if (sepN_length != (time_steps - 1)) {
      time_steps = sepN_length - 1;
      Rf_warningcall(R_NilValue,
        "Resetting time_steps to length of separate_N - 1.");
    }
  }
  
  if (time_steps < 1) throw Rcpp::exception("Option time_steps must be positive.",
    false);
  
  NumericVector output(time_steps + 1);
  double used_phi {0.0};
  double nt {0.0};
  
  output(0) = start_value;
  
  for (int i = 1; i <= time_steps; i++ ) {
    used_phi = 0.0;
    nt = 0.0;
    
    if ((i - time_lag) < 0) {
      if (pre0_subs) {
        used_phi = pre0_value;
      } else {
        used_phi = start_value;
      }
      
      if (base_N) {
        nt = sepN(0);
      } else {
        nt = used_phi;
      }
    } else {
      used_phi = output(i - time_lag);
      
      if (base_N) {
        nt = sepN(i - time_lag);
      } else {
        nt = used_phi;
      }
    }
    
    if (beta > 0.0) {
      if (nt > beta) nt = beta;
    }
    output(i) = used_phi * lambda * (1 - (nt / alpha));
    
    if (substoch > 0) {
      if (output(i) < 0.0) {
        output(i) = 0.0;
      } else if (output(i) > 1.0 && substoch == 2) {
        output(i) = 1.0;
      }
    }
  }
  return output;
}

//' Function to Test Whether a Numeric Vector Is Composed Only of Integers
//' 
//' Function \code{.integer_test()} tests whether the elements of a numeric
//' vector are integers.
//' 
//' @name .integer_test
//' 
//' @param numeric_input A numeric vector to be tested.
//' 
//' @return This function returns an integer equal to the number of elements
//' that are not integers in \code{numeric_input}.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.integer_test)]]
int integer_test(NumericVector numeric_input) {
  int vec_length = numeric_input.length();
  int non_ints = 0;
  
  for (int i = 0; i < vec_length; i++) {
    if (floor(numeric_input(i)) != ceil(numeric_input(i))) {
      non_ints++;
    }
  }
  
  return non_ints;
}

//' Function to Test Whether an Integer Vector Is Entirely Binomial
//' 
//' Function \code{.binomial_test()} tests whether an integer vector is composed
//' only of \code{0} and \code{1} elements.
//' 
//' @name .binomial_test
//' 
//' @param integer_input An integer vector to be tested.
//' 
//' @return This function returns an integer equal to the number of elements
//' that are not \code{0}s and \code{1}s.
//' 
//' @keywords internal
//' @noRd
// [[Rcpp::export(.binomial_test)]]
int binomial_test(NumericVector integer_input) {
  int vec_length = integer_input.length();
  int non_bins = 0;
  
  for (int i = 0; i < vec_length; i++) {
    if (integer_input(i) < 0 || integer_input(i) > 1) {
      non_bins++;
    }
  }
  
  return non_bins;
}

