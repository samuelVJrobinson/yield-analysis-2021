#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density;
  using namespace Eigen; //Needed for sparse structures
  // using namespace atomic; //Needed for logdet
  
  //Simple 1-D gam
  
  //Read data from R------------------
  DATA_VECTOR(Y);               //Response
  DATA_MATRIX(X);               //Design matrix (no intercept)
  DATA_SPARSE_MATRIX(S1);       //Penalty matrix - must be sparse for GMRF to work
  DATA_INTEGER(S1dim);          //Dimensions of S1
  DATA_MATRIX(newX);            //New design matrix for predictions
  
  //Read parameters from R------------
  PARAMETER(b0);       //Intercept
  PARAMETER_VECTOR(smoothCoefs); //Smooth coefficients
  PARAMETER(log_lambda);// Log penalization parameter
  PARAMETER(log_sigma);   // log(SD) of measurement error 
  
  //Transform SD and penalization parameters 
  Type sigma = exp(log_sigma);
  Type lambda = exp(log_lambda);
  
  //Define objective function
  Type nll=0;
  
  // Penalization part ------------------------------------
  
  //logdet requires a matrix type, and GMRF needs a sparse matrix
  // SparseMatrix<Type> S = lambda*S1; //Penalty term * Penalty matrix
  // nll -= 0.5*(logdet(matrix<Type>(S))) - 0.5*(GMRF(S).Quadform(smoothCoefs)); //Returns NaNs. I think this is because logDet(S) is infinite (S = rank deficient?)
  
  //This approach avoids logdet calculation - results are identical as far as I can tell
  nll -= Type(0.5)*S1dim*log_lambda - 0.5*lambda*GMRF(S1).Quadform(smoothCoefs); 
  
  vector<Type> mu = b0 + X*smoothCoefs; //Expected value: mu = Intercept + model matrix %*% coefs
  
  nll -= sum(dnorm(Y, mu, sigma, true)); //Decrement logLik
  // for(int i=0; i<Y.size(); i++){ //Other approach that uses a loop
  //   nll -= dnorm(Y(i), mu(i), sigma, true); 
  // }

  vector<Type> splineForReport = b0 + newX*smoothCoefs; //Generate predictions 
  ADREPORT(splineForReport); //Predicted smoother values
  
  // ADREPORT(b0); //Intercept
  // ADREPORT(smoothCoefs); //Smoothing coefficients - these seem to be left in the model report no matter what
  
  return nll;
}
