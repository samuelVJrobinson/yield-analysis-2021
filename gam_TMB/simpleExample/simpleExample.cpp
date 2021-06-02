#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(NL); //Tooth nitrogen
  DATA_VECTOR(x); //Age
  
  PARAMETER(beta0); //Intercept
  PARAMETER(beta1); //Slope
  PARAMETER(logSigma); //log SD
  
  Type sigma= exp(logSigma); //Get SD
  Type nll = 0; //Set negative log-likelihood to 0
  
  vector<Type> mu = beta0 + x*beta1; //Expected value
  nll -= sum(dnorm(NL,mu,sigma,true)); //Subtract sum of log likelihood from 0

  return(nll);
}
