#include <TMB.hpp>   //Links in the TMB libraries
template<class Type>
Type objective_function<Type>::operator() (){
  using namespace density; //includes GMRF
  using namespace R_inla; //includes SPDE-specific functions, e.g. Q_spde()
  using namespace Eigen; //Needed for sparse structures
                      
  //Read data from R------------------
  DATA_VECTOR(yield); DATA_UPDATE(yield);           // Response
  DATA_VECTOR(logArea); DATA_UPDATE(logArea);       // log-Area of harvested polygons 
  DATA_SCALAR(meanLogArea);                         // Mean of log(area)
  DATA_MATRIX(smoothMat); DATA_UPDATE(smoothMat);   // Design matrix (no intercept)
  DATA_SPARSE_MATRIX(penaltyMat);                   // Penalty matrix - must be sparse for GMRF to work
  DATA_IVECTOR(penaltyDim);                         // Dimensions of penalty sub-matrices
  DATA_MATRIX(predModMat_dist);                     // New smoothing matrices for predictions
  DATA_MATRIX(predModMat_spatial);
  DATA_MATRIX(predModMat_temporal);
  DATA_STRUCT(spatialSPDEMatrices,spde_t);          // Matrices needed for representing the spatial GMRF
  DATA_MATRIX(spatialA); DATA_UPDATE(spatialA);     // Matrix for interpolating points within spatial mesh
  DATA_STRUCT(temporalSPDEMatrices,spde_t);         // Matrices needed for representing the temporal GMRF
  DATA_MATRIX(temporalA); DATA_UPDATE(temporalA);   // Matrix for interpolating points within temporal mesh
  DATA_SCALAR(usePenalty);                          // Should penalty be calculated?

  //Read parameters from R------------
  
  //Mean parameters
  PARAMETER(b0);                              // Intercept
  PARAMETER(b_areaMean);                      // Area (combine speed)
  PARAMETER_VECTOR(smoothCoefsMean);          // Smooth coefficients
  PARAMETER_VECTOR(log_lambdasMean);          // Log penalization parameter
  //Variance parameters
  PARAMETER(b0SD);                            // Intercept
  PARAMETER(b_areaSD);                        // Area (combine speed)
  PARAMETER_VECTOR(smoothCoefsSD);            // Smooth coefficients
  PARAMETER_VECTOR(log_lambdasSD);            // Log penalization parameter
  //SPDE Parameters
  PARAMETER_VECTOR(spdeCoefs_spatial);        // Spatial SPDE coefficients
  PARAMETER(log_tau_spatial);                 // Precision
  PARAMETER(log_kappa_spatial);               // Spatial range
  PARAMETER_VECTOR(spdeCoefs_temporal);       // Spatial SPDE coefficients
  PARAMETER(log_tau_temporal);                // Precision
  PARAMETER(log_kappa_temporal);              // Spatial range
  
  // Transform parameters 
  vector<Type> lambdasMean = exp(log_lambdasMean);
  vector<Type> lambdasSD = exp(log_lambdasSD);
  Type tau_spatial = exp(log_tau_spatial);
  Type kappa_spatial = exp(log_kappa_spatial);
  Type tau_temporal = exp(log_tau_temporal);
  Type kappa_temporal = exp(log_kappa_temporal);
  
  //Recast to sparse matrices
  SparseMatrix<Type> spatialA_sparse = asSparseMatrix(spatialA); 
  SparseMatrix<Type> temporalA_sparse = asSparseMatrix(temporalA);
  
  //Initialize objective function
  Type nll=0;
  
  // Distance, Spatial, and Temporal smoother part ------------------------------------
  
  //Penalizing spline coefficients - taken from pSplines example
  if(usePenalty==1.0){
    int k=0;  // Counter
    for(int i=0;i<penaltyDim.size();i++){ //For each smoother
      int m = penaltyDim(i); //Dimensions of sub-matrix i
      vector<Type> betaMean = smoothCoefsMean.segment(k,m);  //Mean smoother coefficients
      vector<Type> betaSD = smoothCoefsSD.segment(k,m);      //SD smoother coefficients
      
      SparseMatrix<Type> S = penaltyMat.block(k,k,m,m);  // Get penalty sub-matrix for i-th smoother
      
      nll -= 0.5*m*log_lambdasMean(i) - 0.5*lambdasMean(i)*GMRF(S).Quadform(betaMean); //Decrement logLik for each set of coefficients
      nll -= 0.5*m*log_lambdasSD(i) - 0.5*lambdasSD(i)*GMRF(S).Quadform(betaSD);
      
      k += m; //Increment to next sub-matrix
    }
  }
  
  
  //Spatial GMRF part --------------------------------------------------

  //Match location to random intercept, divide by tau (precision parameter) to get spatial random effect
  vector<Type> deltaSpatial = (spatialA_sparse*spdeCoefs_spatial)/tau_spatial;
  //Sparse precision matrix for latent field
  SparseMatrix<Type> Q_spatial = Q_spde(spatialSPDEMatrices,kappa_spatial);
  //NOTE: GMFR returns negative log likelihood, so using += instead of -=
  nll += GMRF(Q_spatial)(spdeCoefs_spatial); //~MVNorm(0, Matern Covariance Matrix) but using SPDE approach

  //Temporal GMRF part -------------------------------------------------

  //Match location to random intercept, divide by tau (precision parameter) to get temporal random effect
  vector<Type> deltaTemporal = (temporalA_sparse*spdeCoefs_temporal)/tau_temporal;
  //Sparse precision matrix for latent field
  SparseMatrix<Type> Q_temporal = Q_spde(temporalSPDEMatrices,kappa_temporal);
  nll += GMRF(Q_temporal)(spdeCoefs_temporal); //~MVNorm(0, Matern Covariance Matrix) but using SPDE approach
  
  // Main model -------------------------------------------

  //Expected value: mu = Intercept + model matrix %*% coefs + spatial/temporal random effects
  vector<Type> mu = b0 + logArea*b_areaMean + smoothMat*smoothCoefsMean + deltaSpatial + deltaTemporal; //Version with SPDE
  vector<Type> logsigma = b0SD + logArea*b_areaSD + smoothMat*smoothCoefsSD; //Version with smoothed SD
  vector<Type> sigma = exp(logsigma);

  for(int i=0; i<yield.size(); i++){
    nll -= dnorm(yield(i), mu(i), sigma(i), true); //Decrement logLik for each data point
  }
  
  // Output for report -------------------------------

  vector<Type> residuals = yield - mu; //Residuals
  for(int i=0; i<yield.size(); i++){
    residuals(i) =  residuals(i)/sigma(i); //Standardize residuals
  }

 //Predictions for smoothers
 // Distance smoothers
 int l = 0; //Counter
 int m = penaltyDim(l); //Which smoother?
 vector<Type> distCoefsMean = smoothCoefsMean.segment(l,m); //Coefficients
 vector<Type> distCoefsSD = smoothCoefsSD.segment(l,m);
 vector<Type> distSmootherMean = b0 + meanLogArea*b_areaMean + predModMat_dist*distCoefsMean; //Generate predictions
 vector<Type> distSmootherSD = b0SD + meanLogArea*b_areaSD + predModMat_dist*distCoefsSD;

  // Spatial smoothers
  l++;
  m = penaltyDim(l);
  vector<Type> spatialCoefsMean = smoothCoefsMean.segment(l,m);
  vector<Type> spatialCoefsSD = smoothCoefsSD.segment(l,m);
  vector<Type> spatialSmootherMean = b0 + meanLogArea*b_areaMean + predModMat_spatial*spatialCoefsMean;
  vector<Type> spatialSmootherSD = b0SD + meanLogArea*b_areaSD + predModMat_spatial*spatialCoefsSD;

  //Temporal smoothers
  l++;
  m = penaltyDim(l);
  vector<Type> temporalCoefsMean = smoothCoefsMean.segment(l,m);
  vector<Type> temporalCoefsSD = smoothCoefsSD.segment(l,m);
  vector<Type> temporalSmootherMean = b0 + meanLogArea*b_areaMean + predModMat_temporal*temporalCoefsMean;
  vector<Type> temporalSmootherSD = b0SD + meanLogArea*b_areaSD + predModMat_temporal*temporalCoefsSD;

  //Distance at which correlation has dropped below 0.1, see p. 4 in Lindgren et al. (2011)
  Type spatialRange = sqrt(8)/kappa_spatial;
  Type temporalRange = sqrt(8)/kappa_temporal;

  REPORT(residuals); //Residuals
  REPORT(spatialRange); //Spatial decay
  REPORT(temporalRange); //Temporal decay
  REPORT(distSmootherMean); //Predicted smoother values
  REPORT(distSmootherSD);
  REPORT(spatialSmootherMean);
  REPORT(spatialSmootherSD);
  REPORT(temporalSmootherMean);
  REPORT(temporalSmootherSD);

  return nll; //Return likelihood
}
