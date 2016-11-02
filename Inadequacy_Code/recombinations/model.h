#ifndef __HYDROGEN_MODEL_H__
#define __HYDROGEN_MODEL_H__

#include "reaction_info.h"
#include <vector>

void
hydrogenComputeModel(
  std::vector<double>& initialValues,
  std::vector<double>& timePoints,
  reaction_info *      rxn,
  std::vector<double>& returnValues);

#endif
