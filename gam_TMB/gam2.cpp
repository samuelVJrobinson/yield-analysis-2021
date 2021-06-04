#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen; //Needed for sparse structures
  
  //1-D gam with 2 predictors
  
  //Read data from R------------------
  DATA_VECTOR(Y);               //Response
  DATA_MATRIX(X);               //Design matrix (no intercept)
  DATA_SPARSE_MATRIX(S1);       //Penalty matrix for smoother 1
  DATA_SPARSE_MATRIX(S2);       //Penalty matrix for smoother 2
  DATA_INTEGER(S1dim);          //Dimensions of S1
  DATA_INTEGER(S2dim);          //Dimensions of S2
  DATA_MATRIX(newX);            //New design matrix for predictions
  
  //Read parameters from R------------
  PARAMETER(b0);       //Intercept
  PARAMETER_VECTOR(smoothCoefs); //Smooth coefficients
  PARAMETER(log_lambda1);// Log penalization parameter
  PARAMETER(log_lambda2);
  PARAMETER(log_sigma);   // log(SD) of measurement error 
  
  
  //Transform SD and penalization parameters 
  Type sigma = exp(log_sigma);
  Type lambda1 = exp(log_lambda1);
  Type lambda2 = exp(log_lambda2);
  
  // Get smoothing coefficients
  vector<Type> smoothCoefs1 = smoothCoefs.segment(0,S1dim);
  vector<Type> smoothCoefs2 = smoothCoefs.segment(S1dim-1,S2dim);
  
  //Define objective function
  Type nll=0;
  
  // Penalization 
  nll -= Type(0.5)*(S1dim*log_lambda1) - 0.5*(lambda1*GMRF(S1).Quadform(smoothCoefs1));
  nll -= Type(0.5)*(S2dim*log_lambda2) - 0.5*(lambda2*GMRF(S2).Quadform(smoothCoefs2)); 

  //Expected value: mu = Intercept + model matrix %*% coefs
  vector<Type> mu = b0 + X*smoothCoefs; 
  
  nll -= sum(dnorm(Y, mu, sigma, true)); //Decrement logLik
  
  // vector<Type> splineForReport = b0 + newX*smoothCoefs; //Generate predictions 
  // ADREPORT(splineForReport); //Predicted smoother values
  
  return nll;
}
