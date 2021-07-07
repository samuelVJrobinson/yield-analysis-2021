#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen; //Needed for sparse structures
  
  //1D GAM with multiple smoothers (using "by" term)
  
  //Read data from R------------------
  DATA_VECTOR(yield);                 //Response
  DATA_VECTOR(logArea);               //log-Area 
  DATA_SCALAR(meanLogArea);           //Mean of log(area)
  DATA_MATRIX(smoothMat);             //Design matrix (no intercept)
  DATA_SPARSE_MATRIX(penaltyMat);     //Penalty matrix - must be sparse for GMRF to work
  DATA_INTEGER(penaltyDim);           //Dimensions of penaltyMat
  DATA_INTEGER(numSmooths);               //Number of smoothers
  DATA_MATRIX(newSmoothMat);          //New smoothing matrix for predictions
  
  //Read parameters from R------------
  PARAMETER(b0);                      //Intercept
  PARAMETER(b_area);                  //Area (combine speed)
  PARAMETER_VECTOR(smoothCoefs);      //Smooth coefficients
  PARAMETER_VECTOR(log_lambdas);      // Log penalization parameter
  PARAMETER(log_sigma);               // log(SD) of measurement error 
  
  //Transform SD and penalization parameters 
  Type sigma = exp(log_sigma);
  vector<Type> lambdas = exp(log_lambdas);
  
  //Initialize objective function
  Type nll=0;
  
  // Penalization part ------------------------------------
  
  for(int i=0;i<numSmooths;i++){ //For each smoother
    //GMRF needs a sparse matrix - This approach avoids logdet calculation
    SparseMatrix<Type> S = lambdas(i)*penaltyMat; //Penalty term * Penalty matrix
    vector<Type> beta_i = smoothCoefs.segment(penaltyDim*i,penaltyDim);       // Coefficients for i-th smoother
    
    nll -= 0.5*penaltyDim*log_lambdas(i) - 0.5*GMRF(S).Quadform(beta_i);  
  }
  
  // Main model -------------------------------------------
  vector<Type> mu = b0 + logArea*b_area + smoothMat*smoothCoefs; //Expected value: mu = Intercept + model matrix %*% coefs
  
  nll -= sum(dnorm(yield, mu, sigma, true)); //Decrement logLik
  
  // Other things to report -------------------------------
  
  vector<Type> splineForReport = b0 + meanLogArea*b_area + newSmoothMat*smoothCoefs; //Generate predictions at mean combine speed
  vector<Type> residuals = yield - mu; //Residuals 
  
  REPORT(residuals); //Residuals
  ADREPORT(splineForReport); //Predicted smoother values
  
  return nll; //Return likelihood
}
