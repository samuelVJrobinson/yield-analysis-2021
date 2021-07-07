// code auto-generated by TMBam
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  using namespace Eigen;
  using namespace R_inla;

  DATA_MATRIX(X); // design matrix
  DATA_VECTOR(y); // response
  // penalty matrices
  DATA_SPARSE_MATRIX(S1);
  DATA_SPARSE_MATRIX(S2);

  PARAMETER(mu); // intercept
  PARAMETER_VECTOR(beta);
  PARAMETER_VECTOR(log_lambda); // smoopar
  PARAMETER(log_sigma);

  // transform smoopars and scale
  vector<Type> lambda = exp(log_lambda);
  Type sigma = exp(log_sigma);

  // split up the random effects coefs
  vector<Type> beta1 = beta.segment(0,9);
  vector<Type> beta2 = beta.segment(9,9);

  // define the real penalties
  Eigen::SparseMatrix<Type> K1 = lambda(0)*S1;
  Eigen::SparseMatrix<Type> K2 = lambda(1)*S2;

  // initialize log-likelihood
  Type nll=0;

  using namespace atomic; // for logdet

  // calculate REML penalty
  nll = -0.5*(logdet(matrix<Type>(K1)) + logdet(matrix<Type>(K2))) +
    0.5*(GMRF(K1).Quadform(beta1) + GMRF(K2).Quadform(beta2));

  // linear predictor
  vector<Type> eta = mu + X*beta;

  for(int i=0; i<y.size(); i++)
      nll -= dnorm(y(i), eta(i), sigma, true);

  return nll;

}