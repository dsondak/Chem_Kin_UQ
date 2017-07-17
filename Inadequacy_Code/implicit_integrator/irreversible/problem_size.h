#ifndef PROBLEM_SIZE_H
#define PROBLEM_SIZE_H

class problem_size
{
  public:
    problem_size (unsigned int n_eq, unsigned int n_species_from_user, unsigned int n_species_d_from_user, 
                  unsigned int n_atoms_from_user, 
                  unsigned int inert_from_user, unsigned int user_n_phis, unsigned int user_n_scenario, 
                  unsigned int user_n_times, unsigned int user_n_reactions, 
                  double user_oxidizer, double user_nitrogen, double user_fuel, bool user_heat_rates, 
                  bool user_init_temperatures, double user_heating_rate, double user_TO, 
                  double user_time_ig, double user_Tig);
    unsigned int n_eq;             // number of equations
    unsigned int n_species;        // number of species (not including inert species)
    unsigned int n_species_d;      // number of species in detailed model (not including inert species)
    unsigned int n_atoms;          // number of atoms in system
    unsigned int n_inert;          // inert species
    unsigned int n_phis;           // number of equivalence ratios
    unsigned int n_scenario;       // number of initial temperatures
    unsigned int n_times;          // number of time points
    unsigned int n_reactions;      // number of reactions in reduced model
    double oxidizer_i;             // Initial concentration of oxidizer (O2)
    double nitrogen;               // Concentration of N2
    double fuel;                   // Stoichiometric factor for fuel (H2)
    bool heat_rates;               // Determine if running heating rate problem
    bool init_temperatures;        // Determine if running initial temperature problem
    double heating_rate;           // Heating rate
    double TO;                     // Initial temperature
    double time_ig;                // Truth data ignition time
    double Tig;                    // Truth data ignition temperature
};

#endif
