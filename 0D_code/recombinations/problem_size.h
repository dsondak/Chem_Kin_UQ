#ifndef PROBLEM_SIZE_H
#define PROBLEM_SIZE_H

#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;

class problem_size
{
  public:
    problem_size (int n_species_from_user, int extra_from_user, int n_reactions_from_user, double heating_rate_from_user, bool include_inad_from_user, int n_atoms_from_user);
    int n_species; // number of species
    int n_extra; // number of extra species (inerts? N2, H2O2)
    int n_reactions; // number of reactions
    double heating_rate; // heating rate
    bool include_inad; // include inadequacy or not
    int n_atoms; // number of additional atoms used in inadequacy operator
    MatrixXd nukj_r; // Reactant stoichiometric coeffs
    MatrixXd nukj_p; // Product stoichiometric coeffs
    MatrixXd nukj; // Total stoichiometric coeffs
    MatrixXi species_to_reaction_map; // Map for species in each reaction
    VectorXd gamma; // Exponent in equilibrium constant
};

#endif
