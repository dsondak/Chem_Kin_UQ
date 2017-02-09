/*-------------------------------------------------------------------
 *
 ***   problem_size.c   ***
 *
 * Class defining the problem size.
 *
 * INPUTS:
 *         n_reactions:  Number of reactions in reduced model
 *
 * MEMBERS:
 *         n_reactions
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

#include "problem_size.h"
#include <iostream>

problem_size::problem_size (unsigned int n_species_from_user, 
                            unsigned int inert_from_user, 
                            unsigned int inad_from_user, 
                            unsigned int equations_from_user, 
                            unsigned int user_n_reactions, 
                            double user_Q)
:
    n_species(n_species_from_user), // set n_species = n_species_from_user
    n_inert(inert_from_user), 
    n_inad(inad_from_user), 
    n_eq(equations_from_user),
    n_reactions(user_n_reactions),
    Q(user_Q), 
    alphas(n_inad, 3), 
    constraints(n_inad)
{
    //alphas(0,0) = -5.87744e+03;
    //alphas(0,1) =  1.75898;
    //alphas(0,2) =  5.48487e-02;
    //alphas(1,0) = -4.66567E+03;
    //alphas(1,1) = 1.38413E+00;
    //alphas(1,2) = 1.20676E-04;

    alphas(0,0) = -5.87744e+03;
    alphas(0,1) =  1.75898;
    alphas(0,2) =  5.48487e-02;
    alphas(1,0) = -5.87744e+03;
    alphas(1,1) =  1.75898;
    alphas(1,2) =  5.48487e-02;

    constraints(0) = 1.63179;
    constraints(1) = 2.72734;
}

