#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla; //includes SPDE-specific functions, e.g. Q_spde()
  using namespace density;
  using namespace Eigen; //Needed for sparse structures
  
  //Read data from R------------------
  DATA_VECTOR(yield);                 //Response
  DATA_STRUCT(spdeMatrices,spde_t);   //Three matrices needed for representing the GMRF, see p. 8 in Lindgren et al. (2011)
  DATA_SPARSE_MATRIX(A);              //Matrix for interpolating points within triangles   
  
  //Read parameters from R------------
  PARAMETER(b0);                  //Intercept
  PARAMETER_VECTOR(spdeCoefs);    //SPDE coefficients (basically random effects coefs)
  PARAMETER(log_tau);             // Tau (Precision parameter)
  PARAMETER(log_kappa);           // Kappa (Spatial decay parameter)
  PARAMETER(log_sigma);           // Sigma (Data variance)
  
  Type tau = exp(log_tau);      // Precision
  Type kappa = exp(log_kappa);  // Spatial decay range
  Type sigma = exp(log_sigma);  // Data variance

  // Spatial interpolation
  vector<Type> delta = (A*spdeCoefs)/tau; //Match location to random intercept, divide by tau (precision parameter) to get spatial random effect
  
  //Construct sparse precision matrix for latent field
  SparseMatrix<Type> Q = Q_spde(spdeMatrices,kappa);

  //Initialize objective function
  Type nll = 0.0;
  
  nll = GMRF(Q)(spdeCoefs); //Basically: ~MVNorm(0, Matern Covariance Matrix) but using SPDE approach
  
  // Main model -------------------------------------------
  vector<Type> mu = b0 + delta; //Expected value: mu = Intercept + spatial random effect
  
  nll -= sum(dnorm(yield, mu, sigma, true)); //Decrement log-likelihood
  
  // Other things to report -------------------------------
  vector<Type> residuals = yield - mu; //Residuals 
  Type range = sqrt(8)/kappa;   //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  
  REPORT(residuals); //Residuals
  ADREPORT(range);
  
  return nll; //Return likelihood
}
