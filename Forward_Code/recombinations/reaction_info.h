#ifndef REACTION_INFO_H
#define REACTION_INFO_H

// User functions
#include "problem_size.h"

//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/reaction_set.h>
#include <antioch/nasa_evaluator.h>
#include <antioch/kinetics_evaluator.h>
//c++
#include <vector>

// define struct that holds all reaction info, except params
struct reaction_info { reaction_info(
  Antioch::ChemicalMixture<double> * chem_mixture,
  Antioch::ReactionSet<double> * reaction_set, 
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * thermo,
  Antioch::KineticsEvaluator<double> * kinetics,
  std::vector<double> & scales,
  int species_from_user,
  int extra_from_user,
  int reactions_from_user,
  double heating_rate_from_user, 
  bool include_inad_from_user, 
  int n_atoms_from_user=0);
 ~reaction_info();

  Antioch::ChemicalMixture<double> * Chem_mixture;
  Antioch::ReactionSet<double> * Reaction_set;
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * Thermo;
  Antioch::KineticsEvaluator<double> * Kinetics;
  std::vector<double> & Scales;
  int n_species;
  int n_extra;
  int n_reactions;
  double heating_rate;
  bool include_inad;
  int n_atoms;
  problem_size ProblemInfo;
};
#endif
