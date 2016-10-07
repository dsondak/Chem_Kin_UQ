/*-------------------------------------------------------------------
 *
 ***   forward_model.c   ***
 *
 * This is the driver to simulate the 0D reactor.
 *
 * INPUTS:
 *         o input.txt (see README on how to construct)
 *         o .xml file containing chemical reactions
 *         o .xml file containing thermodynamic data
 *         o List of species (line 110 of this file) Note that if
 *           the species are changed the code will have to be
 *           recompiled. 
 *
 * OUTPUTS:
 *         o *.h5          contains the time-evolution of each species
 *                         as well as temperature for each scenario (N columns)
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

// User functions
#include "model.h"
#include "reaction_info.h"
#include "write_data.h"
// C libraries
#include "grvy.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <limits>
//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/chemical_species.h>
#include <antioch/chemical_mixture.h>
#include <antioch/reaction_set.h>
#include <antioch/chemkin_parser.h>
#include <antioch/xml_parser.h>
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"
#include <antioch/read_reaction_set_data.h>
// grvy
#include<grvy.h>
#include<sys/time.h>
#include<time.h>

int main()
{

  // Input parameters
  int n_species;           // Number of species (not including extras for N2 and H2O2)
  int n_atoms;             // Number of distinct atoms in the mixture (ignoring N)
  int do_inad;             // Include inadequacy operator
  int n_extra;             // Extra for N2 and H2O2
  int n_phis;              // Equivalence ratios to run
  int n_heating;           // Different heating rates to run
  int n_T;                 // Different starting temperatures to run
  int n_reactions;         // Number of reactions
  int heat_rates;          // Heating rate problem
  int init_temperatures;   // Initial temperature problem
  double timePoint;        // Time-step size
  int numtimes;            // Number of time-steps to run
  double fuel;             // Stoichiometric factor for fuel (H2 here)
  double oxidizer_i;       // Initial concentration of oxidizer
  double nitrogen;         // Concentration of nitrogen
  char *thermo_filename;   // Thermodynamics input file
  char *reaction_filename; // Reaction input file
  char *data_filename;     // Filename to write data to
  
  // Read input file
  grvy_input_fopen("./input.txt"); // Open input file

  grvy_input_fread_int("n_species", &n_species);
  grvy_input_fread_int("n_atoms", &n_atoms);
  grvy_input_fread_int("include_inad", &do_inad);
  grvy_input_fread_int("n_extra", &n_extra);
  grvy_input_fread_int("n_phis", &n_phis);
  grvy_input_fread_int("n_heating", &n_heating);
  grvy_input_fread_int("n_T", &n_T);
  grvy_input_fread_int("n_reactions", &n_reactions);
  grvy_input_fread_int("heat_rates", &heat_rates);
  grvy_input_fread_int("Temperatures", &init_temperatures);
  grvy_input_fread_double("time_points", &timePoint);
  grvy_input_fread_int("num_times", &numtimes);
  grvy_input_fread_double("fuel", &fuel);
  grvy_input_fread_double("oxidizer_i", &oxidizer_i);
  grvy_input_fread_double("nitrogen", &nitrogen);
  grvy_input_fread_char("thermo", &thermo_filename);
  grvy_input_fread_char("reactionset", &reaction_filename);
  grvy_input_fread_char("dataset", &data_filename);

  bool include_inad = false;
  if (do_inad == 1)
  {
     include_inad = true;
  }
  else
  {
     include_inad = false;
  }

  // Scenario parameters
  double phiPoints[n_phis]; // Equivalence ratios
  grvy_input_fread_double_vec("phis", phiPoints, n_phis);

  // Check different possibilities for heating rates and initial temperatures
  // Will either calibrate with initial temperatures or different heating rates
  // but not both.
  if ( (n_heating > 1) && (n_T > 1) )
  {
     std::cout << "You have chosen n_heating > 1 AND n_T > 1.  Only one of these can be > 1." << std::endl;
     exit(0);
  }

  if ( (heat_rates == 1) && (init_temperatures == 1) )
  {
     std::cout << "You cannot run both a heating rate problem and an initial temperature problem.  You must choose one or the other." << std::endl;
     exit(0);
  }

  // Scenario parameters stored as vectors
  double heatPoints[n_heating];
  double tempPoints[n_T];

  grvy_input_fread_double_vec("heating_rate", heatPoints, n_heating);
  grvy_input_fread_double_vec("TO", tempPoints, n_T);

  // Set number of scenario parameters (default = 1)
  int n_scenario = 1;

  // Determine which problem is being run
  bool do_heat_rates = false;
  bool do_initial_temps = false;

  // Modify scenario parameters based on user input
  if (heat_rates == 1)
  {
     n_scenario = n_heating;
     do_heat_rates = true;
  }

  if (init_temperatures == 1)
  {
     n_scenario = n_T;
     do_initial_temps = true;
     heatPoints[0] = 0.0; // No heating rate
  }


  grvy_input_fclose();

  // Total number of variables (species + temperature)
  int dim = n_species + n_extra + 1; // +1 for T
  if (include_inad)
  {
     dim += n_atoms;
  }

  // Initial values of each field
  std::vector<double> initialValues(dim, 0.0);

  // Output file
  //const char fname[] = "truth_data.h5";
  //create_file(fname, numtimes, dim, n_phis*n_scenario);
  create_file(data_filename, numtimes, dim, n_phis*n_scenario);
  int scen; // Current scenario


  /**********************************
  ***
  ***   set up Antioch to do chemistry
  ***
  ***********************************/
       
  // Define species 
  std::vector<std::string> species_str_list;
  species_str_list.reserve(n_species + n_extra);
  species_str_list.push_back("H2");
  species_str_list.push_back("O2");
  species_str_list.push_back("H");
  species_str_list.push_back("O");
  species_str_list.push_back("OH");
  species_str_list.push_back("HO2");
  species_str_list.push_back("H2O");
  species_str_list.push_back("H2O2");
  species_str_list.push_back("N2");

  // Get chemistry for species involved in this reaction
  Antioch::ChemicalMixture<double> chem_mixture( species_str_list );

  // Get thermodynamic data
  Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > nasa_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data( nasa_mixture, thermo_filename, Antioch::XML );

  // Prepare for chemical reactions
  Antioch::ReactionSet    <double> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data_xml<double>( reaction_filename, true, reaction_set );

  // Set up reactions and thermodynamics 
  Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0); // Reactions
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > thermo(nasa_mixture); // Thermodynamics

  std::vector<double> scales(1, 0.); // This is an old variable...can probably remove in the near future

  for (int l = 0; l < n_phis; l++)
  { // Loop over equivalence ratios
      for (int ll = 0; ll < n_scenario; ll++)
      { // Loop over initial temperatures
          // Set initial values of each species and temperature
          initialValues[0] = fuel*phiPoints[l]; // H2
          initialValues[1] = oxidizer_i;        // O2
          initialValues[8] = nitrogen;          // N2
          if (do_initial_temps)
          {
             initialValues[dim-1] = tempPoints[ll]; // initial temperature
          }
          else
          {
             initialValues[dim-1] = tempPoints[0]; // initial temperature
          }

          // Define time
          std::vector<double> timePoints(numtimes,0.);
          for (unsigned int i = 0; i < timePoints.size(); i++)
          {
              timePoints[i] = i * timePoint + 1.0e-06;
          }

          // Initialize solution vector
          std::vector<double> returnValues(dim*numtimes + 2, 0.0); // + 2 for ignition data

          // Call the model
          if (do_heat_rates)
          {
             // This class contains all the chemistry information
             //reaction_info rxnMain(&chem_mixture, &reaction_set, &thermo, &kinetics, scales, n_species, n_extra, heatPoints[ll]);
             reaction_info rxnMain(&chem_mixture, &reaction_set, &thermo, &kinetics, scales, n_species, n_extra, n_reactions, heatPoints[ll], include_inad, n_atoms);
             // Return time points of all species and temperature
             hydrogenComputeModel(initialValues,timePoints,&rxnMain,returnValues);
             // Write out data. 
             scen = l * n_scenario + ll; // Current scenario
             write_file(timePoints, returnValues, data_filename, numtimes, dim, scen, phiPoints[l], heatPoints[ll]);
          }
          else
          {
             // This class contains all the chemistry information
             reaction_info rxnMain(&chem_mixture, &reaction_set, &thermo, &kinetics, scales, n_species, n_extra, n_reactions, heatPoints[0], include_inad, n_atoms);
             // Return time points of all species and temperature
             hydrogenComputeModel(initialValues,timePoints,&rxnMain,returnValues);
             // Write out data. 
             scen = l * n_scenario + ll; // Current scenario
             write_file(timePoints, returnValues, data_filename, numtimes, dim, scen, phiPoints[l], tempPoints[ll]);
          }
 
     } // end scenario loop
  } // end phi loop

  return 0;
}
