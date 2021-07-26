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
  DATA_VECTOR(bStart);               // Start of batches
  DATA_INTEGER(bLength);              // Length of batches
  DATA_INTEGER(nbatches);             // Number of batches
  
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
  int k = 0; //Counter
  int nCoefs = smoothCoefs.size(); //Number of smoothing coefficients
  
  for(int i=0;i<nbatches;i++){ //For each batch
    matrix<Type> sMatBlock = smoothMat.block(k,0,bLength,nCoefs); //Pull out a block of model matrix starting at 0,k going for bLength rows, nCoefs columns
    vector<Type> ySeg = y.segment(k,bLength); //Pull out data for that block
    vector<Type> muSeg = b0 + sMatBlock*smoothCoefs; //Expected value for that block
    
    nll -= sum(dnorm(ySeg, muSeg, sigma, true)); //Decrement logLik
    
    k = k + bLength; //Increment counter
  }
  
  return nll; //Return likelihood
}
