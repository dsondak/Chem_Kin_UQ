/*-------------------------------------------------------------------
 *
 ***   reaction_info.c   ***
 *
 * Structure containing all reaction info, problem_size class and
 * inadequacy_model class.
 *
 * INPUTS:
 *        chem_mixture:  Chemical mixture from Antioch
 *        reaction_set:  Reaction set from Antioch
 *        thermo      :  Thermodynamics data from Antioch
 *        kinetics    :  Reaction kinetics from Antioch
 *        species     :  Number of species in model
 *        inert       :  Additional species
 *        n_reactions :  Number of reactions
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/
#include "reaction_info.h"
#include "problem_size.h"

//Constructor
reaction_info::reaction_info(
    Antioch::ChemicalMixture<double> * chem_mixture,
    Antioch::ReactionSet<double> * reaction_set, 
    Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * thermo,
    Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > * nasa_mixture,
    Antioch::KineticsEvaluator<double> * kinetics,
    std::vector<double> Amat_from_user,
    std::vector<double> bvec_from_user,
    unsigned int species_from_user,
    unsigned int inert_from_user, 
    unsigned int inad_from_user, 
    unsigned int atoms_from_user, 
    unsigned int equations_from_user,
    unsigned int user_n_reactions, 
    double user_Q)
:
  Chem_mixture(chem_mixture),
  Reaction_set(reaction_set),
  Thermo(thermo),
  NASAMixture(nasa_mixture),
  Kinetics(kinetics),
  Amat(Amat_from_user),
  bvec(bvec_from_user),
  n_species(species_from_user),
  n_inert(inert_from_user),
  n_inad(inad_from_user), 
  n_atoms(atoms_from_user), 
  n_eq(equations_from_user),
  n_reactions(user_n_reactions),
  Q(user_Q),
  ProblemInfo(n_species, n_inert, n_inad, n_atoms, n_eq, n_reactions, Q)
{
}

//Destructor
reaction_info::~reaction_info()
{
}
