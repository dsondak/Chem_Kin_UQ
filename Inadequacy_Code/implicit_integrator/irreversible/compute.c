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
#include <grvy.h>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
// user defined functions
#include "compute.h"
#include "likelihood.h"
#include "reaction_info.h"
#include "chemistryVectorRV.h"
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
  int n_atoms;             // Number of distinct atom types in system
  int n_inad;             // Number of atoms
  int n_inert;             // Extra for N2 and H2O2
  int n_species_inad;      // Number of species used in inadequacy model
  int n_species_d;         // Number of species used in detailed model
  int n_phis;              // Equivalence ratios to run
  int n_heating;           // Different heating rates to run
  int n_T;                 // Different starting temperatures to run
  double timePoint;                 // Time-step size
  int n_times;             // Number of time-steps to run
  int n_times_d;           // Number of time-steps in detailed profile
  int n_reactions;         // Number of reactions in mechanism
  int n_reactions_inad;    // Number of reactions in inadequacy model
  double fuel;                      // Stoichiometric factor for fuel (H2 here)
  double oxidizer_i;                // Initial concentration of oxidizer
  double nitrogen;                  // Concentration of nitrogen
  int heat_rates_in;       // Flag for heating rate calibration
  int init_temperature_in; // Flag for initial temperature calibration
  double TO;                        // Initial temperature
  double heating_rate;              // Heating rate
  char *thermo_filename;            // Thermodynamics input file
  char *reaction_filename;          // Reaction input file
  char *data_filename;              // Filename to write data to
  double time_ig;                   // Ignition time from detailed model
  double Tig;                       // Ignition temperature from detailed model

  // Open input parameters file
  grvy_input_fopen("./input.txt");

  // Read in parameters
  grvy_input_fread_int("n_species", &n_species);
  grvy_input_fread_int("n_atoms", &n_atoms);
  grvy_input_fread_int("n_inad", &n_inad);
  grvy_input_fread_int("n_inert", &n_inert);
  grvy_input_fread_int("n_species_inad", &n_species_inad);
  grvy_input_fread_int("n_species_d", &n_species_d);
  grvy_input_fread_int("n_phis", &n_phis);
  grvy_input_fread_int("n_heating", &n_heating);
  grvy_input_fread_int("n_T", &n_T);
  grvy_input_fread_int("heat_rates", &heat_rates_in);
  grvy_input_fread_int("Temperatures", &init_temperature_in);
  grvy_input_fread_double("TO", &TO);
  grvy_input_fread_double("heating_rate", &heating_rate);
  grvy_input_fread_double("time_points", &timePoint);
  grvy_input_fread_int("num_times", &n_times);
  grvy_input_fread_int("n_times_d", &n_times_d);
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

  // Total number of species
  const unsigned int n_species_tot = n_species + n_inert + n_inad;

  // Number of parameters (not hyperparameters)
  const unsigned int n_params = 3 * n_reactions_inad;

  // All of the parameters
  const unsigned int n_params_total = 3 * n_params;

  // Number of equations to solve (species + temperature)
  const unsigned int n_eq = n_species_tot;
 
  // Set number of scenario parameters (default = 1)
  unsigned int n_scenario = 1;

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

  // Read truth data
  const unsigned int n_eq_d = n_species_d + n_inert + 1;
  printf("n_times_d = %2.1i,     n_eq_d = %2.1i\n\n", n_times_d, n_eq_d);
  truth_data detailed_profile("detailed_profile.h5", 1, 1, n_times_d, n_eq_d);

  // Set up data
  std::map<double, double> Tdata;
  for (unsigned int n = 0; n < n_times_d; n++) {
      double time  = detailed_profile.sample_points[n];
      double T     = detailed_profile.observation_data[n* n_eq_d + n_eq_d - 1];
      Tdata[time] = T;
  }


  /*===================================
  ***
  ***   Set up Antioch to do chemistry
  ***
  =====================================*/

  // Define species
  std::vector<std::string> species_str_list;
  species_str_list.reserve(n_species_tot);
  species_str_list.push_back("H2");
  species_str_list.push_back("O2");
  species_str_list.push_back("H");
  species_str_list.push_back("O");
  species_str_list.push_back("OH");
  species_str_list.push_back("HO2");
  species_str_list.push_back("H2O");
  species_str_list.push_back("Hp");
  species_str_list.push_back("Op");
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

  // Number of atoms of each type in the system
  std::vector<double> bvec(n_atoms + 1, 0.0);
  bvec[0] = 2.0 * fuel;
  bvec[1] = 2.0 * oxidizer_i;
  bvec[2] = 2.0 * nitrogen;

  // Species to atoms transition matrix
  std::vector<double> Amat((n_atoms + 1) * n_species_tot, 0.0);
  Amat[0] = 2.0;
  Amat[1] = 0.0;
  Amat[2] = 1.0;
  Amat[3] = 0.0;
  Amat[4] = 1.0;
  Amat[5] = 1.0;
  Amat[6] = 2.0;
  Amat[7] = 1.0;
  Amat[8] = 0.0;
  Amat[9] = 0.0;

  Amat[10] = 0.0;
  Amat[11] = 2.0;
  Amat[12] = 0.0;
  Amat[13] = 1.0;
  Amat[14] = 1.0;
  Amat[15] = 2.0;
  Amat[16] = 1.0;
  Amat[17] = 0.0;
  Amat[18] = 1.0;
  Amat[19] = 0.0;

  Amat[20] = 0.0;
  Amat[21] = 0.0;
  Amat[22] = 0.0;
  Amat[23] = 0.0;
  Amat[24] = 0.0;
  Amat[25] = 0.0;
  Amat[26] = 0.0;
  Amat[27] = 0.0;
  Amat[28] = 0.0;
  Amat[29] = 2.0;

  // Introduce scaling factors for some of the larger parameters
  std::vector<double> scales (n_params, 1.0);

  // Now reset scaling factors for Arrhenius prefactors
  scales[0]  = 1.0e+08;
  scales[3]  = 1.0e+10;
  scales[6]  = 1.0e+08;
  scales[9]  = 1.0e+10;
  scales[12] = 1.0e+07;
  scales[15] = 1.0e+10; 
  scales[18] = 1.0e+07; 
  scales[21] = 1.0e+10; 
  scales[24] = 1.0e+08; 
  scales[27] = 1.0e+10; 

  // Now reset scaling factors for modified Arrhenius parameters
  scales[1]  = 1.0e-02;
  scales[4]  = 1.0;
  scales[7]  = 1.0e-01;
  scales[10] = 1.0;
  scales[13] = 1.0;
  scales[16] = 1.0;
  scales[19] = 1.0e-01;
  scales[22] = 1.0;
  scales[25] = 1.0e-01;
  scales[28] = 1.0;

  // Now reset scaling factors for Arrhenius activation energies
  scales[2]  = 1.0e+05;
  scales[5]  = 1.0e+04;
  scales[8]  = 1.0e+05;
  scales[11] = 1.0e+05;
  scales[14] = 1.0e+05;
  scales[17] = 1.0e+05;
  scales[20] = 1.0e+04;
  scales[23] = 1.0e+05;
  scales[26] = 1.0e+04;
  scales[29] = 1.0e+04;

  // This class contains all the chemistry and problem size information
  reaction_info rxnMain(
      &chem_mixture,
      &reaction_set,
      &thermo,
      &nasa_mixture,
      &kinetics,
      scales,
      Amat, 
      bvec, 
      Tdata, 
      n_eq, 
      n_species,
      n_species_d,
      n_atoms, 
      n_inert,
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

  /*======================================================
  ***
  *** Instantiate the parameter space
  ***
  ========================================================*/
  
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", n_params_total, NULL);

  /*======================================================
  ***
  *** Instantiate the parameter domain
  ***
  ========================================================*/

  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());

  // Domain for all parameters and their hypermeans
  for (int i = 0; i < 2 * n_params; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] =  INFINITY;
  }

  // Domain for hypervariances
  for (int i = 0; i < n_params; i++)
  {
      paramMinValues[2 * n_params + i] = 0.0;
      paramMaxValues[2 * n_params + i] = INFINITY;
  }

  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

  /*======================================================
  ***
  *** Instantiate the likelihood function  object to be 
  *** used by QUESO.
  ***
  ========================================================*/
  
  
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood(env, paramDomain, &rxnMain);

  /*======================================================
  ***
  *** Define the prior
  ***
  ========================================================*/

  // Chemistry prior for all parameters
  QUESO::ChemistryVectorRV<QUESO::GslVector, QUESO::GslMatrix> 
      priorTotal("priorTotal_", paramDomain, rxnMain.InadInfo.n_reactions_inad, rxnMain.InadInfo.n_inad);

  /*======================================================
  ***
  *** Instantiate the inverse problem
  ***
  ========================================================*/

  QUESO::GenericVectorRV<>postTotal("post_", paramSpace);
           
  QUESO::StatisticalInverseProblem<> ip("", NULL, priorTotal, lhood, postTotal); 

  /*======================================================
  ***
  *** Solve the inverse problem
  ***
  ========================================================*/

  // Set up initial parameter values

  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  std::vector<double> 
    init_kinetics = {log(2.2716495127060917), 5.2756035602390422, 1.1021099751123025, 
                     log(3.7714830544483780), 1.0e-16, 8.9433281629646342, 
                     log(1.7721512823487061), 3.8246742586538279, 2.0745500053927372, 
                     log(5.5706916507985680), 1.0e-16, 1.3135622105611983, 
                     log(8.2830919705299377), -1.3585766873971914, 1.3005166005172624, 
                     log(3.5305061634938347), 1.0e-16, 1.7168401051086030,
                     log(6.6857188699084237), 6.6634385904742177, 7.9705662739958876, 
                     log(7.8430783348457687), 1.0e-16, 1.0228335042104413,
                     log(4.3298112708110565), -6.9922795750273214, 9.0747132030759138, 
                     log(3.7484557650026878), 1.0e-16, 4.3359623996979244};

  // First set up the Arrhenius inadequacy parameters
  for (int i = 0; i < n_params; i++)
  {
      paramInitials[i] = init_kinetics[i];
  }

  // Now set up hypermeans of the Arrhenius inadequacy parameters
  for (int i = 0; i < n_params; i++)
  {
      paramInitials[n_params + i] = init_kinetics[i];
  }

  // Finally do the hypervariances
  for (int i = 0; i < n_params; i++)
  {
      paramInitials[2 * n_params + i] = 10.0;
  }

  // Set up covariance matrix for parameters
  QUESO::GslVector diagVec(paramSpace.zeroVector());
  QUESO::GslMatrix proposalCovMatrix(diagVec);

  double var = 1.0e-02;
  for (int i = 0; i < n_params_total; i++)
  {
      if (paramInitials[i] == 1.0e-16)
      {
         proposalCovMatrix(i,i) = var;
      }
      else
      {
         proposalCovMatrix(i,i) = var * paramInitials[i] * paramInitials[i];
      }
      //if (paramInitials[i] != 0.0)
      //{
      //   proposalCovMatrix(i,i) = var * paramInitials[i] * paramInitials[i];
      //}
      //else
      //{
      //   proposalCovMatrix(i,i) = var;
      //}
      //printf("C(%2.1i, %2.1i) = %25.16e     theta = %25.16e\n", i, i, proposalCovMatrix(i,i), paramInitials[i]);
  }
  //exit(0);

  // Solve!
  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);
  
  return;
}
