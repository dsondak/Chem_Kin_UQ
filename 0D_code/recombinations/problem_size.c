// Practice creating the stochastic operator class
#include "problem_size.h"
#include <Eigen/Dense>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;

problem_size::problem_size (int n_species_from_user, int extra_from_user,  
                            int n_reactions_from_user, double heating_rate_from_user, 
                            bool include_inad_from_user, int n_atoms_from_user)
:
    n_species(n_species_from_user), // set n_species = n_species_from_user
    n_extra(extra_from_user), // extra species along for the ride (N2 and H2O2)
    heating_rate(heating_rate_from_user), // heating rate
    include_inad(include_inad_from_user), // include inad. or not
    n_atoms(n_atoms_from_user), // Atoms used in the inadequcy operator
    n_reactions(n_reactions_from_user), // Number of reactions
    nukj_r(n_species + n_atoms, 7), // Reactant stoichiometric coeffs
    nukj_p(n_species + n_atoms, 7),  // Product stoichiometric coeffs
    nukj(n_species + n_atoms, 7), // Total coeffs.
    species_to_reaction_map(7, 3), // Map for species in each reaction
    gamma(7) // Exponent in equilibrium constant
{

    // Reactant stoichiometric coefficients for the catchall reactions
    nukj_r << MatrixXd::Identity(n_species, 7),
              MatrixXd::Zero(n_atoms, 7);

    // Product stoichiometric coefficients for the catchall reactions
    MatrixXd V(n_atoms, 7);

    V << 2, 0, 1, 0, 1, 1, 2,
         0, 2, 0, 1, 1, 2, 1;

    nukj_p << MatrixXd::Zero(n_species, 7), 
              V;

    // Difference between prod. coeffs. and react. coeffs.
    nukj = nukj_p - nukj_r;

    // Exponent on pa/RT, gamma = \sum_{k}{nukj}
    for (int j = 0; j < nukj.cols(); j++)
    {
        gamma(j) = nukj.col(j).sum();
        //std::cout << nukj.col(j).sum() << std::endl;
        //std::cout << gamma(j) << std::endl;
    }

    // Species to reaction map
    int Nc;
    Nc = n_species + n_atoms;
    //MatrixXi species_to_reaction_map(7, 3);
    species_to_reaction_map << 7, 7,  Nc, 
                               8, 8,  Nc,
                               7, Nc, Nc,
                               8, Nc, Nc,
                               7, 8,  Nc,
                               7, 8,  8, 
                               7, 7,  8;

}
