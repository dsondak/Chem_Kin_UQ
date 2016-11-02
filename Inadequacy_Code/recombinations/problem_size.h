#ifndef PROBLEM_SIZE_H
#define PROBLEM_SIZE_H

class problem_size
{
  public:
    problem_size (int n_species_from_user, int n_atoms_from_user, int extra_from_user, 
                  int user_n_phis, int user_n_scenario, int user_n_times, int user_n_reactions, 
                  double user_oxidizer, double user_nitrogen, double user_fuel, bool user_heat_rates, 
                  bool user_init_temperatures, double user_heating_rate, double user_TO, 
                  double user_time_ig, double user_Tig);
    int n_species; // number of species (not including inert species)
    int n_atoms;   // number of atoms in system
    int n_extra;   // inert species
    int n_phis; // number of equivalence ratios
    int n_scenario; // number of initial temperatures
    int n_times; // number of time points
    int n_reactions; // number of reactions in reduced model
    double oxidizer_i; // Initial concentration of oxidizer (O2)
    double nitrogen; // Concentration of N2
    double fuel; // Stoichiometric factor for fuel (H2)
    bool heat_rates; // Determine if running heating rate problem
    bool init_temperatures; // Determine if running initial temperature problem
    double heating_rate;    // Heating rate
    double TO;              // Initial temperature
    double time_ig;         // Truth data ignition time
    double Tig;             // Truth data ignition temperature
};

#endif
