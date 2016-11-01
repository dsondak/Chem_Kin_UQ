/*-------------------------------------------------------------------
 *
 ***   stochastic_operator.c   ***
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

#include "stochastic_operator.h"
#include <iostream>
#include <stdlib.h>
#include <grvy.h>
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;

int gcd(int a, int b)
{
   // This is a function to compute the gcd of integers a and b
   // Taken from:
   // http://codereview.stackexchange.com/questions/66711/greatest-common-divisor
   return b == 0 ? a : gcd(b, a%b);
}

stochastic_operator::stochastic_operator (int n_species_from_user, int n_atoms_from_user, int extra_from_user)
:
    n_species(n_species_from_user), // set n_species = n_species_from_user
    n_atoms(n_atoms_from_user),
    n_extra(extra_from_user),
    species_list(n_species + n_atoms),
    species_list_m(n_species + n_atoms),
    prime_atoms(n_atoms),
    atoms(n_species + n_atoms),
    iatoms(n_species + n_atoms),
    Cmat(n_species + n_atoms, n_species),
    Pmat(n_species, n_species + n_atoms),
    Smat(n_species + n_atoms, n_species + n_atoms)
{

  // Prime number representation of species
  int tmp_vec[n_species + n_atoms]; // use temporary vector b/c grvy doesn't like
                                    // Eigen format or std::vector format
  grvy_input_fopen("./input.txt");
  grvy_input_fread_int_vec("specieslist", tmp_vec, n_species + n_atoms);

  // Now copy temporary vector entries into actual vector
  for (int i = 0; i < n_species + n_atoms; i++)
  {
     species_list(i) = tmp_vec[i];
  }

  // Set up prime number representation of each atom in the system
  int t_vec[n_atoms]; // use temporary vector b/c grvy doesn't like
                      // Eigen format or std::vector format
  grvy_input_fread_int_vec("prime_sieve", t_vec, n_atoms);

  // Now copy temporary vector entries into actual vector
  for (int i = 0; i < n_atoms; i++)
  {
     prime_atoms(i) = t_vec[i];
  }

  grvy_input_fclose();

  // Initialize matrix containing prime factors for each atom and species
  MatrixXd pf = MatrixXd::Zero(n_atoms, n_species + n_atoms);
  int xtest;
  int mod;
  double pexp;

  // species_list_m is the same as species list EXCEPT
  // that it only counts each atom in a species once.
  // So, H20 is represented as 2*2*3 but in
  // species_list_m it is expressed as 2*3.  We call
  // this the zero-multiplicity prime factorization.
  species_list_m.setOnes(n_species+n_atoms);

  for (int k = 0; k < n_species + n_atoms; k++) // Loop over species
  {
     xtest = species_list(k); // Get species number
     for (int i = 0; i < n_atoms; i++) // Loop over atoms
     {
        mod = 0; // Initialize modulus
        while (mod == 0) // While modulus is 0
        {
           mod = xtest % prime_atoms(i); // Get modulus of species and prime sieve
           if (mod == 0)
           {
              // If mod == 0 then we have a prime factor
              xtest = xtest/prime_atoms(i);
              pf(i,k) = pf(i,k) + 1.0; // Increment number of atoms in species
           } // end if
        } // end while mod == 0
        pexp = pf(i,k);
        if (pexp > 1.0)
        {
           pexp = 1.0;
        }
        species_list_m(k) *= pow(prime_atoms(i), pexp);
     } // end loop over atoms
  } // end loop over species

  // Create the C matrix
  double sum; // Needed for column summation of pf
  int gcdij; // Greatest common divisor
  nnz = 0; // Number of nonzero entries in C
  // Set up number of atoms in catchalls (obviously 1)
  for (int i = 0; i < n_atoms; i++)
  {
     atoms(i)  = 1.0;
     iatoms(i) = 1.0;
  }
  double intpart;
  n_eq_nonlin = 0;
  for (int i = 0; i < n_species + n_atoms; i++)
  {
     for (int j = 0; j < n_species; j++)
     {
        if (i < n_atoms) // First n_atoms rows of C matrix
        {
           // Sum column k of pf to find total number
           // of atoms in species
           sum = 0.0;
           for (int k = 0; k < n_atoms; k++)
           {
              sum += pf(k, j + n_atoms);
           }
           atoms(j + n_atoms)  = sum;
           iatoms(j + n_atoms) = 1.0/sum;
           Cmat(i,j) = -pf(i, j + n_atoms)/sum; // Fraction of atom in species
           // Check to see if entry is a fraction (only do for first row)
           // This is to determine number of catchall reactions
           if (i == 0)
           {
              if (modf(Cmat(i,j), &intpart) != 0)
              {
                 n_eq_nonlin += 1; // Increment number of catchall equations
              }
           }
        }
        else if (i == j + n_atoms) // Identity matrix for last rows
        {
           Cmat(i,j) = 1.0;
        }
        else
        {
           Cmat(i,j) = 0.0;
        }
        // Now compute the GCD to determine the number of nonzeros
        gcdij = gcd(species_list_m(j + n_atoms), species_list_m(i));
        if (gcdij >= species_list_m(j + n_atoms))
        {
           nnz += 1;
        }
     } // End loop over n_species
  } // End loop over n_species + n_atoms

  // Create coefficient matrix for catchall reactions
  coeff_catchall.resize(n_atoms, n_eq_nonlin);
  for (int i = 0; i < n_atoms; i++)
  {
      for (int j = 0; j < n_eq_nonlin; j++)
      {
          coeff_catchall(i,j) = pf(i, n_species+n_atoms-n_eq_nonlin+j);
      }
  }

  // Get mapping for nonzero entires
  mapnz.resize(nnz,2);
  int l = 0;
  for (int i = 0; i < n_species; i++)
  {
     for (int j = 0; j < n_species + n_atoms; j++)
     {
        gcdij = gcd(species_list_m(i + n_atoms), species_list_m(j));
        if (gcdij >= species_list_m(i + n_atoms))
        {
           mapnz(l,0) = i;
           mapnz(l,1) = j;
           l += 1;
        }
     }
  }

} // end stochastic_operator constructor

void stochastic_operator::form_operator(VectorXd xi)
{

   MatrixXd Pmat = MatrixXd::Zero(n_species, n_species + n_atoms);

   // Set nonzero values in P matrix
   int i;
   int j;
   for (int k = 0; k < nnz; k++)
   {
      i = mapnz(k,0);
      j = mapnz(k,1);
      Pmat(i,j) = xi[k];
   }

   // Modify some entries in P matrix
   double sum;
   double pre;
   double fact;
   double tmp;
   for (int k = 0; k < nnz; k++)
   {
      // Get nonzero entries
      i = mapnz(k,0);
      j = mapnz(k,1);
      if (j == i + n_atoms) // Then modify entry
      {
         // The next several lines are from Eq. 3.37
         // in Rebecca's dissertation.
         fact = 100.0; // some large value bigger than any expected in C
         for (int n = 0; n < n_atoms; n++)
         // First compute min(abs(Cij)) over rows of col j
         // excluding zero entries
         {
            if (Cmat(n,i) != 0)
               {
                  tmp = std::abs(Cmat(n,i));
                  if (tmp < fact)
                     {
                        fact = tmp;
                     }
               }
         }
         fact = 1.0 / fact;
         sum = 0.0;
         for (int l = 0; l < n_species; l++)
         {
            if (l != i)
            {
               pre = 0.0;
               for (int n = 0; n < n_atoms; n++)
               // Now compute max(abs(Cij)) over rows of col j
               // excluding zero entries
               {
                  if (Cmat(n,l) != 0)
                  {
                     tmp = std::abs(Cmat(n,l));
                     if (tmp > pre)
                     {
                        pre = tmp;
                     }
                  }
               }
               sum += pre * Pmat(l,j);
            }
         }
         // Update the P matrix
         Pmat(i,j) = -Pmat(i,j) - fact * sum;
      } // end flag to modify entry
   } // end loop over nonzero entries

   // Create the stochastic operator matrix
   
   Smat = Cmat * Pmat;

   // std::cout << "S matrix = " << std::endl;
   // std::cout << Smat << std::endl;

}

