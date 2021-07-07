#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace R_inla; //includes SPDE-specific functions, e.g. Q_spde()
  using namespace density;
  using namespace Eigen; //Needed for sparse structures
  
  //Read data from R------------------
  DATA_VECTOR(yield);                        // Response
  DATA_VECTOR(logArea);                      // log-Area 
  DATA_SCALAR(meanLogArea);                  // Mean of log(area)
  DATA_MATRIX(smoothMat);                    // Design matrix (no intercept)
  DATA_SPARSE_MATRIX(penaltyMat);            // Penalty matrix - must be sparse for GMRF to work
  DATA_INTEGER(penaltyDim);                  // Dimensions of penalty matrix
  DATA_MATRIX(newSmoothMat);                 // New smoothing matrix for predictions
  DATA_STRUCT(spatialSPDEMatrices,spde_t);   // Matrices needed for representing the spatial GMRF
  DATA_SPARSE_MATRIX(spatialA);              // Matrix for interpolating points within spatial mesh
  DATA_STRUCT(temporalSPDEMatrices,spde_t);  // Matrices needed for representing the temporal GMRF
  DATA_SPARSE_MATRIX(temporalA);             // Matrix for interpolating points within temporal mesh  
  
  //Read parameters from R------------
  PARAMETER(b0);                              // Intercept
  PARAMETER(b_area);                          //Area (combine speed)
  PARAMETER_VECTOR(smoothCoefs);              // Smooth coefficients
  PARAMETER(log_lambda);                      // Log penalization parameter
  PARAMETER_VECTOR(spdeCoefs_spatial);        // Spatial SPDE coefficients
  PARAMETER(log_tau_spatial);                 // Precision
  PARAMETER(log_kappa_spatial);               // Spatial range
  PARAMETER_VECTOR(spdeCoefs_temporal);       // Spatial SPDE coefficients
  PARAMETER(log_tau_temporal);                // Precision
  PARAMETER(log_kappa_temporal);              // Spatial range
  PARAMETER(log_sigma);                       // Error variance
  
  // Transform parameters 
  Type lambda = exp(log_lambda);
  Type tau_spatial = exp(log_tau_spatial);    
  Type kappa_spatial = exp(log_kappa_spatial);
  Type tau_temporal = exp(log_tau_temporal);    
  Type kappa_temporal = exp(log_kappa_temporal);  
  Type sigma = exp(log_sigma);  
  
  //Initialize objective function
  Type nll=0;
  
  // Distance smoother part ------------------------------------
  
  SparseMatrix<Type> S = lambda*penaltyMat; //Penalty term * Penalty matrix
  // Penalty works with bs = 'ts'/'cs': extra shrinkage
  nll -= 0.5*penaltyDim*log_lambda - 0.5*GMRF(S).Quadform(smoothCoefs);  
  
  //Spatial GMRF part ---------------------------------
  
  //Match location to random intercept, divide by tau (precision parameter) to get spatial random effect
  vector<Type> deltaSpatial = (spatialA*spdeCoefs_spatial)/tau_spatial; 
  //Sparse precision matrix for latent field
  SparseMatrix<Type> Q_spatial = Q_spde(spatialSPDEMatrices,kappa_spatial);   
  //NOTE: GMFR returns negative log likelihood, so using += instead of -=
  nll += GMRF(Q_spatial)(spdeCoefs_spatial); //~MVNorm(0, Matern Covariance Matrix) but using SPDE approach
  
  //Temporal GMRF part ---------------------------------
  
  //Match location to random intercept, divide by tau (precision parameter) to get temporal random effect
  vector<Type> deltaTemporal = (temporalA*spdeCoefs_temporal)/tau_temporal; 
  //Sparse precision matrix for latent field
  SparseMatrix<Type> Q_temporal = Q_spde(temporalSPDEMatrices,kappa_temporal);   
  
  nll += GMRF(Q_temporal)(spdeCoefs_temporal); //~MVNorm(0, Matern Covariance Matrix) but using SPDE approach
  
  // Main model -------------------------------------------
  
  //Expected value: mu = Intercept + model matrix %*% coefs + spatial/temporal random effects
  vector<Type> mu = b0 + logArea*b_area + smoothMat*smoothCoefs + deltaSpatial + deltaTemporal;
  
  // nll -= sum(dnorm(yield, mu, sigma, true)); //Decrement negative log likelihood
  
  vector<Type> residuals = yield - mu; //Residuals 
  for(int i=0; i<yield.size(); i++){
    nll -= dnorm(yield(i), mu(i), sigma, true); //Decrement logLik for each data point
    residuals(i) =  (yield(i) - mu(i))/sigma; //Calculate residuals
  }
  
  // Things to report -------------------------------
  
  vector<Type> splineForReport = b0 + meanLogArea*b_area + newSmoothMat*smoothCoefs; //Generate predictions
  
  //Distance at which correlation has dropped to 0.1, see p. 4 in Lindgren et al. (2011)
  Type spatialRange = sqrt(8)/kappa_spatial;   
  Type temporalRange = sqrt(8)/kappa_temporal;   
  
  REPORT(residuals); //Residuals
  ADREPORT(spatialRange);
  ADREPORT(temporalRange);
  ADREPORT(splineForReport); //Predicted smoother values
  
  return nll; //Return likelihood
}
