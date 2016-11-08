/*-------------------------------------------------------------------
 *
 ***   compute.c   ***
 *
 * This is the driver to perform model inadequacy on the stochastic
 * operator.
 *
 * This file does two things:
 *      1.)  Perform the statistical inverse problem to calibrate
 *           uncertain parameters.  These include the reaction rate
 *           parameters and the stochastic operator parameters.
 *
 *      2.)  Perform the statistical forward problem to get the QOI.
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

// C libraries
#include <stdio.h>
#include <stdlib.h>
#include <grvy.h>
// user defined functions
#include "compute.h"
#include "likelihood.h"
#include "reaction_info.h"
#include "model.h"
#include "write_data.h"
//#include "chemistryVectorRV.h"
//antioch
#include <antioch/vector_utils.h>
#include <antioch/antioch_asserts.h>
#include <antioch/chemical_species.h>
#include <antioch/chemical_mixture.h>
#include <antioch/cea_mixture.h>
#include <antioch/cea_evaluator.h>
#include <antioch/reaction_set.h>
#include <antioch/xml_parser.h>
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_mixture_parsing.h"
#include <antioch/read_reaction_set_data.h>
//queso
#include <queso/GslVector.h>
#include <queso/GslMatrix.h>
#include <queso/GenericScalarFunction.h>
#include <queso/GenericVectorFunction.h>
#include <queso/GenericVectorRV.h>
#include <queso/GaussianVectorRV.h>
#include <queso/LogNormalVectorRV.h>
#include <queso/UniformVectorRV.h>
#include <queso/StatisticalInverseProblem.h>
#include <queso/StatisticalForwardProblem.h>
#include <sys/time.h>
#include <cmath>
// Eigen functions
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using Eigen::MatrixXd;
using Eigen::VectorXd;

void computeAllParams(const QUESO::FullEnvironment& env) {


  /*==================================
  ***
  ***  Preamble output 
  ***
  ====================================*/
  
  struct timeval timevalNow;
  
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "\n Executing SIP and SFP "
              << ctime(&timevalNow.tv_sec)
              << "\n my fullRank = "         << env.fullRank()
              << "\n my subEnvironmentId = " << env.subId()
              << "\n my subRank = "          << env.subRank()
              << "\n my interRank = "        << env.inter0Rank()
               << std::endl << std::endl;
  }

  // More printing 
  if ((env.subDisplayFile()) && (env.displayVerbosity() >= 2))
  {
     *env.subDisplayFile() << "Run starting at: " << ctime(&timevalNow.tv_sec) << std::endl;
  }
  env.fullComm().Barrier();
  env.subComm().Barrier();
  
  /*================================================================
  ***
  ***   Statistical inverse problem (SIP)
  ***   Calibrate parameters. 
  ***
  ================================================================*/
  gettimeofday(&timevalNow, NULL);
  if (env.fullRank() == 0) {
    std::cout << "Beginning 'SIP -> all parameters estimation' at "
              << ctime(&timevalNow.tv_sec)
              << std::endl;
  }

  // Input parameters
  int n_species;           // Number of species (not including extras for N2 and H2O2)
  int n_atoms;             // Number of atoms
  int n_extra;             // Extra for N2 and H2O2
  int n_species_inad;      // Number of species used in inadequacy model
  int n_phis;              // Equivalence ratios to run
  int n_heating;           // Different heating rates to run
  int n_T;                 // Different starting temperatures to run
  double timePoint;        // Time-step size
  int n_times;             // Number of time-steps to run
  int n_reactions;         // Number of reactions in mechanism
  int n_reactions_inad;    // Number of reactions in inadequacy model
  double fuel;             // Stoichiometric factor for fuel (H2 here)
  double oxidizer_i;       // Initial concentration of oxidizer
  double nitrogen;         // Concentration of nitrogen
  int heat_rates_in;       // Flag for heating rate calibration
  int init_temperature_in; // Flag for initial temperature calibration
  double TO;               // Initial temperature
  double heating_rate;     // Heating rate
  char *thermo_filename;   // Thermodynamics input file
  char *reaction_filename; // Reaction input file
  char *data_filename;     // Filename to write data to
  double time_ig;          // Ignition time from detailed model
  double Tig;              // Ignition temperature from detailed model

  // Open input parameters file
  grvy_input_fopen("./input.txt");

  // Read in parameters
  grvy_input_fread_int("n_species", &n_species);
  grvy_input_fread_int("n_atoms", &n_atoms);
  grvy_input_fread_int("n_extra", &n_extra);
  grvy_input_fread_int("n_species_inad", &n_species_inad);
  grvy_input_fread_int("n_phis", &n_phis);
  grvy_input_fread_int("n_heating", &n_heating);
  grvy_input_fread_int("n_T", &n_T);
  grvy_input_fread_int("heat_rates", &heat_rates_in);
  grvy_input_fread_int("Temperatures", &init_temperature_in);
  grvy_input_fread_double("TO", &TO);
  grvy_input_fread_double("heating_rate", &heating_rate);
  grvy_input_fread_double("time_points", &timePoint);
  grvy_input_fread_int("num_times", &n_times);
  grvy_input_fread_int("num_reactions", &n_reactions);
  grvy_input_fread_int("n_reactions_inad", &n_reactions_inad);
  grvy_input_fread_double("fuel", &fuel);
  grvy_input_fread_double("oxidizer_i", &oxidizer_i);
  grvy_input_fread_double("nitrogen", &nitrogen);
  grvy_input_fread_char("thermo", &thermo_filename);
  grvy_input_fread_char("reactionset", &reaction_filename);
  grvy_input_fread_double("time_ig", &time_ig);
  grvy_input_fread_char("dataset", &data_filename);
  grvy_input_fread_double("Tig", &Tig);

  // Close input parameters file
  grvy_input_fclose();

  // Check different possibilities for heating rates and initial temperatures
  // Will either calibrate with initial temperatures or different heating rates
  // but not both.
  if ( (n_heating > 1) && (n_T > 1) )
  {
     std::cout << "You have chosen n_heating > 1 AND n_T > 1.  Only one of these can be > 1." << std::endl;
     exit(0);
  }

  // Set number of scenario parameters (default = 1)
  int n_scenario = 1;

  // Modify scenario parameters based on user input
  if (n_heating > 1)
  {
     n_scenario = n_heating;
  }

  if (n_T > 1)
  {
     n_scenario = n_T;
  }

  bool heat_rates = false;
  bool init_temperatures = false;

  if ((heat_rates_in == 1) && (init_temperature_in == 1))
  {
     std::cout << "Calibration over heating rates and initial temperatures simultaneously is prohibited." << std::endl;
     exit(0);
  }

  if (heat_rates_in == 1)
  {
     heat_rates = true;
  }

  if (init_temperature_in == 1)
  {
     init_temperatures = true;
     heating_rate = 0.0;
  }

  // Number of parameters in reaction terms
  // Usually will be three (A, b, Ea) but may
  // need to modify if including three-body
  // reactions.
  const unsigned int n_ks = 3 * n_reactions;

  // Total number of variables (species + temperature)
  const int dim = n_species + n_atoms + n_extra + 1; // + 1 for T

  /*===================================
  ***
  ***   Set up Antioch to do chemistry
  ***
  =====================================*/

  // Read in chemical reactions

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
  species_str_list.push_back("N2");

  // Get chemistry for species involved in this reaction
  Antioch::ChemicalMixture<double> chem_mixture( species_str_list );

  // Get thermodynamic data
  Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > nasa_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data( nasa_mixture, thermo_filename, Antioch::XML );

  // Prepare for chemical reactions
  Antioch::ReactionSet<double> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data_xml<double>(reaction_filename, true, reaction_set );

  // Set up reaction and thermodynamics 
  Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0);
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > thermo(nasa_mixture);

  // This class contains all the chemistry and problem size information
  reaction_info rxnMain(
      &chem_mixture,
      &reaction_set,
      &thermo,
      &kinetics,
      n_species,
      n_atoms,
      n_extra,
      n_species_inad,
      n_phis,
      n_scenario,
      n_times,
      n_reactions,
      n_reactions_inad,
      oxidizer_i,
      nitrogen,
      fuel,
      heat_rates,
      init_temperatures,
      heating_rate,
      TO,
      time_ig,
      Tig);

  const unsigned int Np = 3;    // parameters, hypermean, hyperstd
  const unsigned int n_arr = 3; // number of Arrhenius parameters
  const unsigned int order = 2; // polynomial order of catchall thermochemistry

  const unsigned int n_params = Np * n_arr * rxnMain.ProblemInfo.n_reactions + 
                                Np * rxnMain.ProblemInfo.n_atoms + 
                                Np * rxnMain.ProblemInfo.n_atoms * order + 
                                Np * rxnMain.ProblemInfo.n_atoms;

/*  RUN THE MODEL */

/*

  create_file(data_filename, n_times, dim, n_phis*n_scenario);

  rxnMain.inad_model.alphas(0,0) = -3.90372558e+04;
  rxnMain.inad_model.alphas(0,1) =  13.6559654;
  rxnMain.inad_model.alphas(0,2) =  1.20459536e-03;

  rxnMain.inad_model.alphas(1,0) = -3.90372558e+04;
  rxnMain.inad_model.alphas(1,1) =  13.6559654;
  rxnMain.inad_model.alphas(1,2) =  1.20459536e-03;

  rxnMain.inad_model.betas(0) = -17.0857829376;
  rxnMain.inad_model.betas(1) = -17.0857829376;

  rxnMain.inad_model.Aj(0) = 1.0e+09;
  rxnMain.inad_model.Aj(1) = 2.0e+05;
  rxnMain.inad_model.Aj(2) = 1.0e+09;
  rxnMain.inad_model.Aj(3) = 1.5e+03;
  rxnMain.inad_model.Aj(4) = 1.17e+03;

  rxnMain.inad_model.bj(0) = 0.5;
  rxnMain.inad_model.bj(1) = 0.75;
  rxnMain.inad_model.bj(2) = 0.5;
  rxnMain.inad_model.bj(3) = 0.25;
  rxnMain.inad_model.bj(4) = 0.1;

  rxnMain.inad_model.Ej(0) = 1.65e+05;
  rxnMain.inad_model.Ej(1) = 1.65e+05;
  rxnMain.inad_model.Ej(2) = 1.65e+05;
  rxnMain.inad_model.Ej(3) = 1.65e+05;
  rxnMain.inad_model.Ej(4) = 1.65e+05;

  std::vector<double> initial_conditions(dim, 0.0);

  initial_conditions[0] = 1.0 * rxnMain.ProblemInfo.fuel;
  initial_conditions[1] = rxnMain.ProblemInfo.oxidizer_i;
  initial_conditions[7] = rxnMain.ProblemInfo.nitrogen;
  initial_conditions[dim-1] = rxnMain.ProblemInfo.TO;

  std::vector<double> sample_points(n_times, 0.0);
  for (int n = 0; n < n_times; n++)
  {
      sample_points[n] = n * timePoint + 1.0e-06;
  }

  std::vector<double> returnValues(n_times * dim + 2, 0.0); // Added two for ignition data

  hydrogenComputeModel(initial_conditions,sample_points,&rxnMain,returnValues);

  int scen = 0;
  write_file(sample_points, returnValues, data_filename, n_times, dim, scen, 1.0, heating_rate);
*/

/**** CHECKING INADEQUACY CLASS. ****

  double temperature = 2000.0;
  
  rxnMain.inad_model.thermo(temperature);

  std::vector<double> Yinad(rxnMain.inad_model.n_species_inad, 0.0);

  Yinad[0] = 2.0;
  Yinad[1] = 1.0;

  double R_universal = Antioch::Constants::R_universal<double>();
  double RT = R_universal * temperature;

  int Nc = rxnMain.ProblemInfo.n_species;

  Antioch::TempCache<double> temp_cache(temperature);

  // Set up s/R - h_RT for each species
  std::vector<double> delta_k(Nc, 0.0);

  // First just do H2 and O2
  for (int k = 0; k < 2; k++)
  {   
      delta_k[k] = rxnMain.Thermo->s_over_R(temp_cache, k) - 
                   rxnMain.Thermo->h_over_RT(temp_cache, k);
  }
  // Skipped H and O so need to decrement index by 2
  for (int k = 4; k < rxnMain.ProblemInfo.n_species; k++)
  {   
      delta_k[k - 2] = rxnMain.Thermo->s_over_R(temp_cache, k) - 
                   rxnMain.Thermo->h_over_RT(temp_cache, k);
  }
  // Now do catchalls (still need to decrement by 2)
  for (int k = 0; k < n_atoms; k++)
  {
      delta_k[rxnMain.ProblemInfo.n_species + k - 2] = rxnMain.inad_model.s_prime(k) / R_universal - 
                                                   rxnMain.inad_model.h_prime(k) / RT;
  }

  rxnMain.inad_model.progress_rate(Yinad, temperature, R_universal, delta_k);

  std::cout << rxnMain.inad_model.rj << std::endl;
*/

  return;
}
