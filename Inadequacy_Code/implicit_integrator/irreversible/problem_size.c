/*-------------------------------------------------------------------
 *
 ***   problem_size.c   ***
 *
 * Class defining the problem size.
 *
 * INPUTS:
 *         n_phis:       Number of equivalence ratios
 *         n_scenario:         Number of initial temperatures
 *         n_reactions:  Number of reactions in reduced model
 *         oxidizer:     Initial concentration of oxidizer
 *         nitrogen:     Initial concentration of nitrogen
 *         fuel:         Initial concentration of fuel
 *
 * MEMBERS:
 *         n_phis
 *         n_scenario
 *         n_times
 *         n_reactions
 *         oxidizer
 *         nitrogen
 *         fuel
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

#include "problem_size.h"
#include <iostream>

problem_size::problem_size (unsigned int n_eq_from_user, unsigned int n_species_from_user, 
                            unsigned int n_species_d_from_user, 
                            unsigned int n_atoms_from_user, unsigned int inert_from_user, 
                            unsigned int user_n_phis, unsigned int user_n_scenario, 
                            unsigned int user_n_times, unsigned int user_n_reactions,
                            double user_oxidizer, double user_nitrogen, double user_fuel, 
                            bool user_heat_rates, bool user_init_temperatures, 
                            double user_heating_rate, double user_TO, double user_time_ig, 
                            double user_Tig)
:
    n_eq(n_eq_from_user),
    n_species(n_species_from_user),
    n_species_d(n_species_d_from_user),
    n_atoms(n_atoms_from_user),
    n_inert(inert_from_user),
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
    Tig(user_Tig)
{
}
