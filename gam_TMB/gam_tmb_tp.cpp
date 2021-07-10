#include <TMB.hpp>  
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen; 
  
  //2D thin-plate regression spline from mgcv
  
  //Read data from R------------------
  DATA_VECTOR(y);                     // Response
  DATA_MATRIX(smoothMat);             // Design matrix (no intercept)
  DATA_SPARSE_MATRIX(penaltyMat);     // Penalty matrix
  DATA_INTEGER(penaltyDim);           // Dimensions of penalty matrix
  
  //Read parameters from R------------
  PARAMETER(b0);                      // Intercept
  PARAMETER_VECTOR(smoothCoefs);      // Smooth coefficients
  PARAMETER(log_lambda);              // Log penalization parameter
  PARAMETER(log_sigma);               // Log(SD) of measurement error 
  
  //Transform SD and penalization parameters 
  Type sigma = exp(log_sigma);
  Type lambda = exp(log_lambda);
  
  //Initialize objective function
  Type nll=0;
  
  // Penalization ------------------------------------
  
  SparseMatrix<Type> S = lambda*penaltyMat; //Penalty term * Penalty matrix
  nll -= 0.5*penaltyDim*log_lambda - 0.5*GMRF(S).Quadform(smoothCoefs);
  
  // Main model -------------------------------------------
  vector<Type> mu = b0 + smoothMat*smoothCoefs; //Expected value
  
  nll -= sum(dnorm(y, mu, sigma, true)); //Decrement logLik
  
  return nll; //Return likelihood
}
