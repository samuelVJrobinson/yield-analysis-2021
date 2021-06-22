#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen; //Needed for sparse structures
  
  //Simple 1-D gam
  
  //Read data from R------------------
  DATA_VECTOR(yield);                 //Response
  DATA_MATRIX(smoothMat);             //Design matrix (no intercept)
  DATA_SPARSE_MATRIX(penaltyMat);     //Penalty matrix - must be sparse for GMRF to work
  DATA_INTEGER(penaltyDim);           //Dimensions of penaltyMat
  DATA_MATRIX(newSmoothMat);          //New smoothing matrix for predictions

  //Read parameters from R------------
  PARAMETER(b0);       //Intercept
  PARAMETER_VECTOR(smoothCoefs); //Smooth coefficients
  PARAMETER(log_lambda);// Log penalization parameter
  PARAMETER(log_sigma);   // log(SD) of measurement error 
  
  //Transform SD and penalization parameters 
  Type sigma = exp(log_sigma);
  Type lambda = exp(log_lambda);
  
  //Initialize objective function
  Type nll=0;
  
  // Penalization part ------------------------------------
  
  //GMRF needs a sparse matrix - This approach avoids logdet calculation
  SparseMatrix<Type> S = lambda*penaltyMat; //Penalty term * Penalty matrix
  // Penalty that works with bs = 'ts'/'cs': extra shrinkage
  nll -= 0.5*penaltyDim*log_lambda - 0.5*GMRF(S).Quadform(smoothCoefs);  

  // Main model -------------------------------------------
  vector<Type> mu = b0 + smoothMat*smoothCoefs; //Expected value: mu = Intercept + model matrix %*% coefs
  
  nll -= sum(dnorm(yield, mu, sigma, true)); //Decrement logLik
  
  
  // Other things to report -------------------------------
  
  vector<Type> splineForReport = b0 + newSmoothMat*smoothCoefs; //Generate predictions
  vector<Type> residuals = yield - mu; //Residuals 
  
  REPORT(residuals); //Residuals
  ADREPORT(splineForReport); //Predicted smoother values
  
  return nll; //Return likelihood
}
