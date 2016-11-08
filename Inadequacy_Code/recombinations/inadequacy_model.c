/*-------------------------------------------------------------------
 *
 ***   inadequacy_model.c   ***
 *
 * Stochastic operator class.
 *
 * INPUTS:
 *         n_species:  Number of species in model
 *         n_atoms:    Number of atoms in model
 *         n_extra:    Number of additional species in model
 *
 * MEMBERS:
 *         n_species
 *         n_atoms
 *         n_extra
 *         n_eq_nonlin   :  Number of equations for catchall reactions
 *         nnz           :  Number of nonzeros in stochastic operator
 *         coeff_catchall:  Coefficient matrix for catchall reactions
 *         Smat          :  Created in method form_operator (line 205)
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

#include "inadequacy_model.h"
#include "reaction_info.h"
#include <iostream>
#include <stdlib.h>
#include <grvy.h>
// Eigen functions
#include <Eigen/Dense>
//antioch
#include <antioch/kinetics_evaluator.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

inadequacy_model::inadequacy_model (int n_species_from_user, int n_atoms_from_user, int extra_from_user, 
                                    int n_species_inad_from_user, int n_reactions_inad_from_user)
:
    n_species(n_species_from_user),               // set n_species = n_species_from_user
    n_atoms(n_atoms_from_user),                   // number of atoms in the system
    n_extra(extra_from_user),                     // number of inert species
    n_species_inad(n_species_inad_from_user),     // number of inadequacy species
    n_reactions_inad(n_reactions_inad_from_user), // number of inadequacy reactions
    nukj_r(n_species_inad, n_reactions_inad),     // Reactant stoichiometric coeffs
    nukj_p(n_species_inad, n_reactions_inad),     // Product stoichiometric coeffs
    nukj(n_species_inad, n_reactions_inad),       // Total coeffs.
    gamma(n_reactions_inad),                      // Exponent in equilibrium constant
    alphas(n_atoms, 3),                           // Enthalpy coefficients for virtual species
    betas(n_atoms),                               // Entropy coefficients for virtual species
    Aj(n_reactions_inad),                         // Arrhenius prefactor for inadequacy model
    bj(n_reactions_inad),                         // Modified Arrhenius exponent for inad. model
    Ej(n_reactions_inad),                         // Activation energy for inadequacy model
    h_prime(n_atoms),                             // Enthalpy for virtual species
    cp_prime(n_atoms),                            // Specific heat for virtual species
    s_prime(n_atoms),                             // Entropy for virtual species
    rj(n_reactions_inad)                          // Progress rate for inadequacy model
{

    // Reactant stoichiometric coefficients for the catchall reactions
    nukj_r << 1, 0, 0, 0, 0, 
              0, 1, 0, 0, 0,
              0, 0, 1, 0, 0, 
              0, 0, 0, 1, 0, 
              0, 0, 0, 0, 1, 
              0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0;

    // Product stoichiometric coefficients for the catchall reactions
    nukj_p << 0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 
              0, 0, 0, 0, 0, 
              2, 0, 1, 1, 2, 
              0, 2, 1, 2, 1;

    // Difference between prod. coeffs. and react. coeffs.
    nukj = nukj_p - nukj_r;

    // Exponent on pa/RT, gamma = \sum_{k}{nukj}
    for (int j = 0; j < nukj.cols(); j++)
    {
        gamma(j) = nukj.col(j).sum();
    }

} // end inadequacy_model constructor

//void inadequacy_model::thermo(MatrixXd alphas, VectorXd betas, double T)
void inadequacy_model::thermo(double T)
{

    for (unsigned int m = 0; m < n_atoms; m++)
    {
        h_prime(m)  = alphas(m,0) + alphas(m,1) * T + alphas(m,2) * T * T;
        cp_prime(m) =               alphas(m,1)  + 2.0 * alphas(m,2) * T;
        s_prime(m)  = betas(m) + alphas(m,1) * log(T) + 2.0 * alphas(m,2) * T;
    }

}

void inadequacy_model::progress_rate(std::vector<double> Yinad, double T, double R, 
                                     std::vector<double> delta_k)
{

     Antioch::TempCache<double> temp_cache(T);

     double RT = R * T;
     double pa = 1.0e+05; // 1 bar in Pa
     double pa_RT = pa / RT;
     double exp_arg_j;

/*
     // Make up some Arrhenius coeffs. for now
     VectorXd A(n_reactions_inad);
     VectorXd beta(n_reactions_inad);
     VectorXd Ea(n_reactions_inad);

     A    << 1.0e+09, 2.0e+05, 1.0e+09, 1.5e+03, 1.17e+03;
     beta << 0.5, 0.75, 0.5, 0.25, 0.1;
     Ea   << 1.65e+05, 1.65e+05, 1.65e+05, 1.65e+05, 1.65e+05;
*/
   
     double kfj; // forward reaction rate coeff.
     double kej; // equilibrium constant
     double kbj; // backward reaction rate coeff.
     double rfj; // forward progress rate
     double rbj; // backward progress rate

     for (int j = 0; j < n_reactions_inad; j++)
     {
         kfj = Aj(j) * pow(T, bj(j)) * exp(-Ej(j) / RT);
         exp_arg_j = 0.0;
         for (int k = 0; k < n_species_inad; k++)
         {
             exp_arg_j  += nukj(k,j) * delta_k[k];
         }
         kej = pow(pa_RT, gamma(j)) * exp(exp_arg_j);
         kbj = kfj / kej;
         rfj = 1.0;
         rbj = 1.0;
         // Calculate reaction rates
         for (int k = 0; k < n_species_inad; k++)
         {
             rfj *= pow(Yinad[k], nukj_r(k,j));
             rbj *= pow(Yinad[k], nukj_p(k,j));
         }
         rj(j) = kfj * rfj - kbj * rbj;
     }

}

