#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla; //includes SPDE-specific functions, e.g. Q_spde()
  using namespace density;
  using namespace Eigen; //Needed for sparse structures
  
  // //Earlier version using AR1 process. Fits, but impossible to get reports due to large sample size
  // //Read data from R------------------
  // DATA_VECTOR(yield);                 //Response
  // 
  // //Read parameters from R------------
  // PARAMETER(b0);                  // Intercept
  // PARAMETER_VECTOR(u);            // Latent AR1 process
  // PARAMETER(atanh_phi);           // AR1 parameter
  // PARAMETER(log_sigmaU);          // Variance for AR1 process
  // PARAMETER(log_sigma);           // Variance for error
  // 
  // //Transform parameters
  // Type phi = tanh(atanh_phi);      
  // Type sigmaU = exp(log_sigmaU);   
  // Type sigma = exp(log_sigma);   
  // 
  // //Initialize objective function
  // Type nll = 0.0;
  // 
  // //AR1 process x sigmaU = latent temporal random effect
  // nll = SCALE(AR1(phi), sigmaU)(u); 
  // 
  // // Main model -------------------------------------------
  // vector<Type> mu = b0 + u; //Expected value: mu = Intercept + temporal random effect
  // 
  // nll -= sum(dnorm(yield, mu, sigma, true)); //Decrement log-likelihood
  // 
  // // Other things to report -------------------------------
  // vector<Type> residuals = yield - mu; //Residuals 
  // REPORT(residuals); //Residuals
  // ADREPORT(u);
  
  return nll; //Return likelihood
}
