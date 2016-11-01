#ifndef __HYDROGEN_LIKELIHOOD_H__
#define __HYDROGEN_LIKELIHOOD_H__

#include "reaction_info.h"
#include "truth_data.h"
#include <queso/GslMatrix.h>
#include <queso/ScalarFunction.h>

template<class V, class M>
class Likelihood : public QUESO::BaseScalarFunction<V, M>
{
public:

  Likelihood(const QUESO::BaseEnvironment& env, const QUESO::VectorSet<V, M> & domainSet, reaction_info * rxnInfo);

  const QUESO::BaseEnvironment* m_env;
  reaction_info * m_rxnMain;
  int n_phis;
  int n_times;
  int n_scen;
  int num_fields;
  int n_species;
  truth_data obs_data;

  virtual double lnValue(
         const QUESO::GslVector& paramValues,
         const QUESO::GslVector* paramDirection,
         QUESO::GslVector*       gradVector,
         QUESO::GslMatrix*       hessianMatrix,
         QUESO::GslVector*       hessianEffect) const;
  
  virtual double actualValue(
         const QUESO::GslVector& paramValues,
         const QUESO::GslVector* paramDirection,
         QUESO::GslVector*       gradVector,
         QUESO::GslMatrix*       hessianMatrix,
         QUESO::GslVector*       hessianEffect) const;

  ~Likelihood();
};

#endif
