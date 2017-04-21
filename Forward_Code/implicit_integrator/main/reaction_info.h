#ifndef REACTION_INFO_H
#define REACTION_INFO_H

// User functions
#include "problem_size.h"
//c++
#include <vector>
//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/reaction_set.h>
#include <antioch/nasa_evaluator.h>
#include <antioch/nasa_mixture.h>
#include <antioch/kinetics_evaluator.h>

// define struct that holds all reaction info, except params
struct reaction_info {
  reaction_info(
      Antioch::ChemicalMixture<double> * chem_mixture,
      Antioch::ReactionSet<double> * reaction_set, 
      Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * thermo,
      Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > * nasa_mixture,
      Antioch::KineticsEvaluator<double> * kinetics,
      std::vector<double> Amat_from_user,
      std::vector<double> bvec_from_user,
      unsigned int species_from_user,
      unsigned int inert_from_user, 
      unsigned int inad_from_user, 
      unsigned int atoms_from_user, 
      unsigned int equations_from_user,
      unsigned int user_n_reactions, 
      double user_Q);
 ~reaction_info();

  Antioch::ChemicalMixture<double> * Chem_mixture;
  Antioch::ReactionSet<double> * Reaction_set;
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * Thermo;
  Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > * NASAMixture;
  Antioch::KineticsEvaluator<double> * Kinetics;
  std::vector<double> Amat;
  std::vector<double> bvec;
  unsigned int n_species;
  unsigned int n_inert;
  unsigned int n_inad;
  unsigned int n_atoms;
  unsigned int n_eq;
  unsigned int n_reactions;
  double Q;
  problem_size ProblemInfo;
};
#endif
