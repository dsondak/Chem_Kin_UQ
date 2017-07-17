/*-------------------------------------------------------------------
 *
 ***   inadequacy_model.c   ***
 *
 * Stochastic operator class.
 *
 * INPUTS:
 *         n_species:  Number of species in model
 *         n_inad:    Number of atoms in model
 *         n_inert:    Number of additional species in model
 *
 * MEMBERS:
 *         n_species
 *         n_inad
 *         n_inert
 *         n_eq_nonlin   :  Number of equations for catchall reactions
 *         nnz           :  Number of nonzeros in stochastic operator
 *         coeff_catchall:  Coefficient matrix for catchall reactions
 *         Smat          :  Created in method form_operator (line 205)
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

#include "inadequacy_model.h"

inadequacy_model::inadequacy_model (unsigned int n_inad_from_user, 
                                    unsigned int n_reactions_inad_from_user)
:
    n_inad(n_inad_from_user),                    // number of atoms in the system
    n_reactions_inad(n_reactions_inad_from_user) // number of inadequacy reactions
{
}

