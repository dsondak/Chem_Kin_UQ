#ifndef __HYDROGEN_QOI_H__
#define __HYDROGEN_QOI_H__

#include "reaction_info.h"
#include <queso/GslMatrix.h>
#include <queso/DistArray.h>

struct qoiRoutine_Data
{
  qoiRoutine_Data(const QUESO::BaseEnvironment& env,
      const std::vector<double> & phis,
      const std::vector<double> & scenario,
      const std::vector<double> & times,
      std::vector<double> & concs,
      reaction_info * rxnInfo);
 ~qoiRoutine_Data();

  const QUESO::BaseEnvironment* m_env;
  const std::vector<double> & m_phis;   //
  const std::vector<double> & m_scenario;   //
  const std::vector<double> & m_times;   //
  const std::vector<double> & m_concs;   // sjfj
  reaction_info       * m_rxnMain;
};

void qoiRoutine(
  const QUESO::GslVector&                     paramValues,
  const QUESO::GslVector*                     paramDirection,
  const void*                                 functionDataPtr,
        QUESO::GslVector&                     qoiValues,
        QUESO::DistArray<QUESO::GslVector* >* gradVectors,
        QUESO::DistArray<QUESO::GslMatrix* >* hessianMatrices,
        QUESO::DistArray<QUESO::GslVector* >* hessianEffects);

#endif
