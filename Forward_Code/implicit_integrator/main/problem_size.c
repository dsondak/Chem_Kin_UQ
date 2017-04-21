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
                            unsigned int atoms_from_user, 
                            unsigned int equations_from_user, 
                            unsigned int user_n_reactions, 
                            double user_Q)
:
    n_species(n_species_from_user), // set n_species = n_species_from_user
    n_inert(inert_from_user), 
    n_inad(inad_from_user), 
    n_atoms(atoms_from_user), 
    n_eq(equations_from_user),
    n_reactions(user_n_reactions),
    Q(user_Q)
{
}

