/*
 * Compute reaction rates using the reduced 
 * mechanism and data from the detailed and 
 * reduced model solutions.
 *
 * Also compute the Jacobian of the reaction 
 * rates at each time.
 *
 * Save everything to HDF5 files.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>

// user defined functions
#include "write_data.h"
#include "truth_data.h"

/* Header files with a description of contents used */

//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/chemical_species.h>
#include <antioch/chemical_mixture.h>
#include <antioch/reaction_set.h>
#include <antioch/chemkin_parser.h>
#include <antioch/xml_parser.h>
#include <antioch/kinetics_evaluator.h>
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"
#include <antioch/read_reaction_set_data.h>

/*********************************
 *
 * Main Program
 *
 *********************************/

int main()
{

  bool detailed = false;
  bool reduced  = true;

  unsigned int Nt = 50000; // # of time-steps
  unsigned int Ns = 7;     // Reduced species
  unsigned int n_eq = 9;   // Make reduced default

  // Total number of fields in detailed and reduced models
  // Includes species, N2 and temperature
  if (detailed) {
     n_eq = 10;
  }
  else if (reduced) {
     n_eq = 9;
  }
  else {
     printf("Model not recognized:  Assuming reduced.");
  }

  /***********************************
  * Read in data and compute RHS 
  * with Antioch using detailed and 
  * reduced data sets.
  ***********************************/

  // Read in data from HDF5 files
  truth_data solution("reduced_solution.h5", 1, 1, Nt, n_eq);

  // Moles of N2 is constant so just set it here
  double x_N2 = solution.observation_data[n_eq - 2];

  double heating_rate = 5.0e+06; // Heating rate

  /***********************************
   * Set up reactions with Antioch
   ***********************************/

  // Define species 
  std::vector<std::string> species_str_list;
  species_str_list.reserve(Ns + 1); // +1 for N2
  species_str_list.push_back("H2");
  species_str_list.push_back("O2");
  species_str_list.push_back("H");
  species_str_list.push_back("O");
  species_str_list.push_back("OH");
  species_str_list.push_back("HO2");
  species_str_list.push_back("H2O");
  species_str_list.push_back("N2");

  // Get chemistry for species involved in this reaction
  Antioch::ChemicalMixture<double> chem_mixture( species_str_list );

  //// Thermodynamics
  Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > nasa_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data( nasa_mixture, "nasa7_thermo_reduced.xml", Antioch::XML );
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > thermo(nasa_mixture); // Thermodynamics

  // Reactions
  Antioch::ReactionSet<double> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data_xml<double>( "reduced_mech_5.xml", true, reaction_set );
  Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0); // Reactions

  // Set some constants
  double R_universal = Antioch::Constants::R_universal<double>();
  double press = 1.0e+05;

  std::vector<double> xs(Ns, 0.0); // Species (in moles)
  double T; // Temperature solution

  std::vector<double> fx(Ns*Nt, 0.0);
  std::vector<double> fT(Nt, 0.0);

  // Loop over each time
  for (unsigned int n = 0; n < Nt; n++) {
      // Begin by reading in solution
      for (unsigned int k = 0; k < Ns; k++) {
          xs[k] = solution.observation_data[n_eq*n + k];
      }
      // Get temperature
      T = solution.observation_data[n_eq*n + n_eq - 1];

      /* Compute volume from ideal gas law */

      // First get the total moles in the system
      double total_moles = 0.0;
      for (unsigned int k = 0; k < Ns; k++) {
          total_moles += xs[k];
      }

      // Add in other species to get the total moles
      if (detailed) {
         double x_H2O2 = solution.observation_data[n_eq*n + n_eq - 3];
         total_moles += x_H2O2 + x_N2;
      }
      else if (reduced) {
         total_moles += x_N2;
      }

      // Finally compute volume
      Antioch::TempCache<double> temp_cache(T);
      double V = total_moles * R_universal * T / press;

      // Molar concentrations
      std::vector<double> molar_densities(Ns+1,0);
      for (unsigned int k = 0; k < Ns; k++) {
          molar_densities[k] = xs[k] / V;
      }
      molar_densities[Ns] = x_N2 / V;

      // Get some thermodynamic data for each species
      std::vector<double> h_RT_minus_s_R(Ns+1, 0.0);
      thermo.h_RT_minus_s_R(temp_cache, h_RT_minus_s_R);

      // Derivative of Gibbs free energy for each species wrt T
      std::vector<double> dh_RT_minus_s_R_dT(Ns+1, 0.0);
      thermo.dh_RT_minus_s_R_dT(temp_cache, dh_RT_minus_s_R_dT);

      // Derivatives of RHS wrt T and species
      std::vector<double> dmole_sources_dT(Ns+1, 0.0);
      std::vector<std::vector<double> > dmole_sources_dX_s(Ns+1);

      for (unsigned int k = 0; k < Ns+1; k++) {
          dmole_sources_dX_s[k].resize(Ns+1);
      }

      // Perform kinetics calculations 
      std::vector<double> mole_sources(Ns+1, 0.0);
      kinetics.compute_mole_sources_and_derivs(
          T,                  // temperature needed for calculations
          molar_densities,    // current concentrations (in moles)
          h_RT_minus_s_R,     // exponent in equilibrium constant for reverse rates
          dh_RT_minus_s_R_dT, // derivative of Keq. exp.
          mole_sources,       // RHS
          dmole_sources_dT,   // Derivative of RHS wrt temperature
          dmole_sources_dX_s  // Derivative of RHS wrt species
      );

      // RHS for species molar concentrations
      // (Not saving RHS for N2)
      for (unsigned int k = 0; k < Ns; k++) {
          fx[Ns*n + k] = mole_sources[k] * V;
      }

      // Energy equation
      double h_dot_mole_sources = 0.0;
      double cp_dot_species     = 0.0;
      double molecular_weight   = 0.0;

      for (unsigned int k = 0; k < Ns; k++) {
          molecular_weight    = R_universal / chem_mixture.R(k);
          h_dot_mole_sources += thermo.h(temp_cache, k) * molecular_weight * mole_sources[k];
          cp_dot_species     += thermo.cp(temp_cache, k) * molecular_weight * xs[k];
      }

      fT[n] = (-h_dot_mole_sources * V + heating_rate) / cp_dot_species;

  } // end time loop

  create_file("reduced_rhs.h5", Nt, Ns);
  write_file(fx, fT, "reduced_rhs.h5", Nt, Ns);

  return(0);

}

