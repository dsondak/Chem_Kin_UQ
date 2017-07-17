#ifndef REACTION_INFO_H
#define REACTION_INFO_H

// User functions
#include "problem_size.h"
#include "inadequacy_model.h"
#include "LinearInterpolation.h"
//c++
#include <vector>
//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/reaction_set.h>
#include <antioch/nasa_evaluator.h>
#include <antioch/nasa_mixture.h>
#include <antioch/nasa_mixture_parsing.h>
#include <antioch/kinetics_evaluator.h>

// define struct that holds all reaction info, except params
struct reaction_info {
  reaction_info(
      Antioch::ChemicalMixture<double> * chem_mixture,
      Antioch::ReactionSet<double> * reaction_set, 
      Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * thermo,
      Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > * nasa_mixture,
      Antioch::KineticsEvaluator<double> * kinetics,
      std::vector<double> & scales,
      std::vector<double> Amat_from_user,
      std::vector<double> bvec_from_user,
      std::map<double, double> Tdata_from_user, 
      unsigned int n_eq_from_user,
      unsigned int species_from_user,
      unsigned int species_d_from_user,
      unsigned int atoms_from_user,
      unsigned int inert_from_user,
      unsigned int species_inad_from_user,
      unsigned int user_n_phis,
      unsigned int user_n_scenario,
      unsigned int user_n_times,
      unsigned int user_n_reactions,
      unsigned int reactions_inad_from_user,
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
  Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > * NASAMixture;
  Antioch::KineticsEvaluator<double> * Kinetics;
  std::vector<double> & scale_factors;
  std::vector<double> Amat;
  std::vector<double> bvec;
  std::map<double, double> Tdata;
  unsigned int n_eq;
  unsigned int n_species;
  unsigned int n_species_d;
  unsigned int n_atoms;
  unsigned int n_inert;
  unsigned int n_species_inad;
  unsigned int n_phis;
  unsigned int n_scenario;
  unsigned int n_times;
  unsigned int n_reactions;
  unsigned int n_reactions_inad;
  double oxidizer_i;
  double nitrogen;
  double fuel;
  bool heat_rates;
  bool init_temperatures;
  double heating_rate;
  double TO;
  double time_ig;
  double Tig;
  problem_size ProblemInfo;
  inadequacy_model InadInfo;
  LinearInterpolation LinInterp;
};
#endif
