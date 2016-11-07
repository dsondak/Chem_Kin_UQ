#ifndef INADEQUACY_MODEL_H
#define INADEQUACY_MODEL_H

#include <vector>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class inadequacy_model
{
  public:
    inadequacy_model (int n_species_from_user, int n_atoms_from_user, int extra_from_user, 
                      int n_species_inad_from_user, int n_reactions_inad_from_user);
    int n_species;     // number of species
    int n_atoms;       // number of atoms in the system
    int n_extra;       // number of extra species (inerts? N2, H2O2)
    int n_species_inad; // number of inadequacy species
    int n_reactions_inad; // number of inadequacy reactions
    MatrixXd nukj_r;   // reactant stoich. coeffs.
    MatrixXd nukj_p;   // product stoich. coeffs.
    MatrixXd nukj;     // nukj_p - nukj_r
    VectorXd gamma;    // exponent in equilibrium constant
    VectorXd h_prime;  // Enthalpy for virtual species
    VectorXd cp_prime; // Specific heat for virtual species
    VectorXd s_prime;  // Entropy for virtual species
    VectorXd rj;       // Progress rate
    void thermo(MatrixXd alphas, VectorXd betas, double T);
    void progress_rate(std::vector<double> Yinad, double T, double R, 
                       std::vector<double> delta_k);
};

#endif
