/*-------------------------------------------------------------------
 *
 ***   reaction_info.c   ***
 *
 * Structure containing all reaction info, problem_size class and
 * classes.
 *
 * INPUTS:
 *        chem_mixture:  Chemical mixture from Antioch
 *        reaction_set:  Reaction set from Antioch
 *        thermo      :  Thermodynamics data from Antioch
 *        kinetics    :  Reaction kinetics from Antioch
 *        species     :  Number of species in model
 *        atoms       :  Number of atoms in model
 *        inert       :  Additional species
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
#include "problem_size.h"
#include "inadequacy_model.h"

//Constructor
reaction_info::reaction_info(
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
    double user_Tig)
:
  Chem_mixture(chem_mixture),
  Reaction_set(reaction_set),
  Thermo(thermo),
  NASAMixture(nasa_mixture),
  Kinetics(kinetics),
  scale_factors(scales),
  Amat(Amat_from_user),
  bvec(bvec_from_user),
  Tdata(Tdata_from_user), 
  n_eq(n_eq_from_user), 
  n_species(species_from_user),
  n_species_d(species_d_from_user),
  n_atoms(atoms_from_user),
  n_inert(inert_from_user),
  n_species_inad(species_inad_from_user),
  n_phis(user_n_phis),
  n_scenario(user_n_scenario),
  n_times(user_n_times),
  n_reactions(user_n_reactions),
  n_reactions_inad(reactions_inad_from_user),
  oxidizer_i(user_oxidizer),
  nitrogen(user_nitrogen),
  fuel(user_fuel),
  heat_rates(user_heat_rates),
  init_temperatures(user_init_temperatures),
  heating_rate(user_heating_rate),
  TO(user_TO),
  time_ig(user_time_ig),
  Tig(user_Tig),
  InadInfo(n_species_inad, n_reactions_inad),
  ProblemInfo(n_eq, n_species, n_species_d, n_atoms, 
              n_inert, n_phis, n_scenario, n_times, 
              n_reactions, oxidizer_i, nitrogen, 
              fuel, heat_rates, init_temperatures, 
              heating_rate, TO, time_ig, Tig),
  LinInterp(Tdata)
{
}

//Destructor
reaction_info::~reaction_info()
{
}
