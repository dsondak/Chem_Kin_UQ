#ifndef __KINETICS_SOLVE_H__
#define __KINETICS_SOLVE_H__

#include "reaction_info.h"
#include <vector>

void
kinetics_forward(
  std::vector<double>& initialValues,
  std::vector<double>& timePoints,
  reaction_info *      rxn,
  std::vector<double>& returnValues);

#endif
