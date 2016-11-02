#ifndef REACTION_INFO_H
#define REACTION_INFO_H

// User functions
#include "inadequacy_model.h"
#include "problem_size.h"
//c++
#include <vector>
//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/reaction_set.h>
#include <antioch/nasa_evaluator.h>
#include <antioch/kinetics_evaluator.h>

// define struct that holds all reaction info, except params
struct reaction_info {
  reaction_info(
      Antioch::ChemicalMixture<double> * chem_mixture,
      Antioch::ReactionSet<double> * reaction_set, 
      Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * thermo,
      Antioch::KineticsEvaluator<double> * kinetics,
      int species_from_user,
      int atoms_from_user,
      int extra_from_user,
      int user_n_phis,
      int user_n_scenario,
      int user_n_times,
      int user_n_reactions,
      double user_oxidizer,
      double user_nitrogen,
      double user_fuel,
      bool user_heat_rates,
      bool user_init_temperatures,
      double user_heating_rate,
      double user_TO,
      double user_time_ig,
      double user_Tig);
 ~reaction_info();

  Antioch::ChemicalMixture<double> * Chem_mixture;
  Antioch::ReactionSet<double> * Reaction_set;
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * Thermo;
  Antioch::KineticsEvaluator<double> * Kinetics;
  int n_species;
  int n_atoms;
  int n_extra;
  int n_phis;
  int n_scenario;
  int n_times;
  int n_reactions;
  double oxidizer_i;
  double nitrogen;
  double fuel;
  bool heat_rates;
  bool init_temperatures;
  double heating_rate;
  double TO;
  double time_ig;
  double Tig;
  std::vector<double> model_params;
  inadequacy_model inad_model;
  problem_size ProblemInfo;
};
#endif
