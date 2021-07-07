#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(obs);
  PARAMETER(phi);
  PARAMETER(logSigma);
  PARAMETER_VECTOR(u);
  using namespace density;
  Type sigma= exp(logSigma);
  Type f = 0;
  f += SCALE(AR1(phi), sigma)(u);
  vector<Type> p = u.exp()/(1.0+u.exp()); //invlogit(u);
  f -= dbinom(obs, Type(1), p, true).sum();
  ADREPORT(sigma);
  ADREPORT(p);
  return f;
  
  //Example from https://www.kaggle.com/berent/template-model-builder-tmb
}
