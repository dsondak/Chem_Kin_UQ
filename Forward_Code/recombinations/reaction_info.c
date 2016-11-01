#include "reaction_info.h"
#include "problem_size.h"

//Constructor
reaction_info::reaction_info(
    Antioch::ChemicalMixture<double> * chem_mixture,
    Antioch::ReactionSet<double> * reaction_set, 
    Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > * thermo,
    Antioch::KineticsEvaluator<double> * kinetics,
    std::vector<double> & scales,
    int species_from_user,
    int extra_from_user,
    int reactions_from_user,
    double heating_rate_from_user, 
    bool include_inad_from_user, 
    int n_atoms_from_user)
: 
  Chem_mixture(chem_mixture),
  Reaction_set(reaction_set),
  Thermo(thermo),
  Kinetics(kinetics),
  Scales(scales),
  n_species(species_from_user),
  n_extra(extra_from_user),
  n_reactions(reactions_from_user),
  heating_rate(heating_rate_from_user),
  include_inad(include_inad_from_user),
  n_atoms(n_atoms_from_user),
  ProblemInfo(n_species, n_extra, n_reactions, heating_rate, include_inad, n_atoms)
{
}

//Destructor
reaction_info::~reaction_info()
{
}
