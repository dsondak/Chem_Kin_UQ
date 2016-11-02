/*-------------------------------------------------------------------
 *
 ***   reaction_info.c   ***
 *
 * Structure containing all reaction info, problem_size class and
 * inadequacy_model class.
 *
 * INPUTS:
 *        chem_mixture:  Chemical mixture from Antioch
 *        reaction_set:  Reaction set from Antioch
 *        thermo      :  Thermodynamics data from Antioch
 *        kinetics    :  Reaction kinetics from Antioch
 *        species     :  Number of species in model
 *        atoms       :  Number of atoms in model
 *        extra       :  Additional species
 *        n_phis      :  Number of equivalence ratios
 *        n_scenario        :  Number of initial temperatures
 *        n_times     :  Number of sample times
 *        n_reactions :  Number of reactions
 *        oxidizer    :  Initial concentration of oxidizer
 *        nitrogen    :  Initial concentration of nitrogen
 *        fuel        :  Initial concentration of fuel
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/
#include "reaction_info.h"
#include "inadequacy_model.h"
#include "problem_size.h"

//Constructor
reaction_info::reaction_info(
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
    double user_Tig)
:
  Chem_mixture(chem_mixture),
  Reaction_set(reaction_set),
  Thermo(thermo),
  Kinetics(kinetics),
  n_species(species_from_user),
  n_atoms(atoms_from_user),
  n_extra(extra_from_user),
  n_phis(user_n_phis),
  n_scenario(user_n_scenario),
  n_times(user_n_times),
  n_reactions(user_n_reactions),
  oxidizer_i(user_oxidizer),
  nitrogen(user_nitrogen),
  fuel(user_fuel),
  heat_rates(user_heat_rates),
  init_temperatures(user_init_temperatures),
  heating_rate(user_heating_rate),
  TO(user_TO),
  time_ig(user_time_ig),
  Tig(user_Tig),
  model_params(),
  inad_model(n_species, n_atoms, n_extra),
  ProblemInfo(n_species, n_atoms, n_extra, 
              n_phis, n_scenario, n_times, 
              n_reactions, oxidizer_i, nitrogen, 
              fuel, heat_rates, init_temperatures, 
              heating_rate, TO, time_ig, Tig)
{
}

//Destructor
reaction_info::~reaction_info()
{
}
