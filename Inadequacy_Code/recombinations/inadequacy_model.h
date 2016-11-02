#ifndef STOCHASTIC_OPERATOR_H
#define STOCHASTIC_OPERATOR_H

#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
using Eigen::VectorXi;

class stochastic_operator
{
  public:
    stochastic_operator (int n_species_from_user, int n_atoms_from_user, int extra_from_user);
    int n_atoms; // number of atoms in the system
    int n_species; // number of species
    int n_extra; // number of extra species (inerts? N2, H2O2)
    VectorXi species_list; // prime number representation of each species
    VectorXi species_list_m; // 0-multiplicity prime number repr of each species (probably redundant)
    VectorXi prime_atoms; // Prime number presentation of each atom in the system
    MatrixXd Cmat;
    MatrixXd Pmat;
    MatrixXd Smat;
    MatrixXd coeff_catchall;
    MatrixXi mapnz;
    VectorXd atoms;
    VectorXd iatoms;
    int nnz; // number of nonzeros in C matrix
    int n_eq_nonlin; // number of catchall reactions
    void form_operator(VectorXd);
};

#endif
