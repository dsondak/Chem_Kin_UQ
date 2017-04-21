#ifndef PROBLEM_SIZE_H
#define PROBLEM_SIZE_H

class problem_size
{
  public:
    problem_size (unsigned int n_species_from_user, 
                  unsigned int inert_from_user,
                  unsigned int inad_from_user,
                  unsigned int atoms_from_user,
                  unsigned int equations_from_user, 
                  unsigned int user_n_reactions, 
                  double user_Q);
    unsigned int n_species;   // number of species (not including inert species)
    unsigned int n_inert;     // inert species
    unsigned int n_inad;      // inadequacy species
    unsigned int n_atoms;     // number of atoms types
    unsigned int n_eq;        // number of equations to solve
    unsigned int n_reactions; // number of reactions in reduced model
    double       Q;           // heating rate
};

#endif
