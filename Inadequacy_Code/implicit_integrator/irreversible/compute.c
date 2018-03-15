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
//#include <grvy.h>
#include <cmath>
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

// YAML
#include <yaml-cpp/yaml.h>

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

  int my_rank = env.fullRank(); 

  // Load file
  std::string input_file = "input_p" + std::to_string(my_rank) + ".yaml";
  //std::cout << "input file name is:  " << input_file << std::endl << std::endl;
  YAML::Node input = YAML::LoadFile(input_file);

  // read data
  int n_species = input["n_species"].as<int>();
  printf("n_species = %2.1i\n", n_species);
  int n_atoms = input["n_atoms"].as<int>();
  printf("n_atoms = %2.1i\n", n_atoms);
  int n_inad = input["n_inad"].as<int>();
  printf("n_inad = %2.1i\n", n_inad);
  int n_inert = input["n_inert"].as<int>();
  printf("n_inert = %2.1i\n", n_inert);
  int n_species_inad = input["n_species_inad"].as<int>();
  printf("n_species_inad = %2.1i\n", n_species_inad);
  int n_species_d = input["n_species_d"].as<int>();
  printf("n_species_d = %2.1i\n", n_species_d);
  int n_phis = input["n_phis"].as<int>();
  printf("n_phis = %2.1i\n", n_phis);
  int n_heating = input["n_heating"].as<int>();
  printf("n_heating = %2.1i\n", n_heating);
  int n_T = input["n_T"].as<int>();
  printf("n_T = %2.1i\n", n_T);
  double timePoint = input["time_points"].as<double>();
  printf("timePoint = %25.16f\n", timePoint);
  int n_times = input["num_times"].as<int>();
  printf("n_times = %2.1i\n", n_times);
  int n_times_d = input["n_times_d"].as<int>();
  printf("n_times_d = %2.1i\n", n_times_d);
  int n_reactions = input["num_reactions"].as<int>();
  printf("n_reactions = %2.1i\n", n_reactions);
  int n_reactions_inad = input["n_reactions_inad"].as<int>();
  printf("n_reactions_inad = %2.1i\n", n_reactions_inad);
  double fuel = input["fuel"].as<double>();
  printf("fuel = %25.16f\n", fuel);
  double oxidizer_i = input["oxidizer_i"].as<double>();
  printf("oxidizer_i = %25.16f\n", oxidizer_i);
  double nitrogen = input["nitrogen"].as<double>();
  printf("nitrogen = %25.16f\n", nitrogen);
  int heat_rates_in = input["heat_rates"].as<int>();
  printf("heat_rates_in = %2.1i\n", heat_rates_in);
  int init_temperature_in = input["Temperatures"].as<int>();
  printf("init_temperature_in = %2.1i\n", init_temperature_in);
  double TO = input["TO"].as<double>();
  printf("TO = %25.16f\n", TO);
  double heating_rate = input["heating_rate"].as<double>();
  printf("heating_rate = %25.16f\n", heating_rate);
  double time_ig = input["time_ig"].as<double>();
  printf("time_ig = %25.16f\n", time_ig);
  double Tig = input["Tig"].as<double>();
  printf("Tig = %25.16f\n", Tig);

  std::string thermo_filename = input["thermo_filename"].as<std::string>();
  printf("thermo_filename = %s\n", thermo_filename.c_str());
  std::string reaction_filename = input["reaction_filename"].as<std::string>();
  printf("reaction_filename = %s\n", reaction_filename.c_str());
  std::string data_filename = input["data_filename"].as<std::string>();
  printf("data_filename = %s\n", data_filename.c_str());

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

  // Get the file name for the exact solution
  char file_name[23];
  char proc_id[2];
  strcpy(file_name, "detailed_profile_p");
  sprintf(proc_id, "%d", my_rank);
  strcat(file_name, proc_id);
  strcat(file_name, ".h5");

  // Get file name for observation data
  char data_fname[17];
  strcpy(data_fname, "truth_data_p");
  sprintf(proc_id, "%d", my_rank);
  strcat(data_fname, proc_id);
  strcat(data_fname, ".h5");

  printf("%s\n", data_fname);

  truth_data detailed_profile(file_name, 1, 1, n_times_d, n_eq_d);

  // Set up data
  std::map<double, double> Tdata;
  for (int n = 0; n < n_times_d; n++) {
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
  for (unsigned int i = 0; i < 2 * n_params; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] =  INFINITY;
  }

  // Domain for hypervariances
  for (unsigned int i = 0; i < n_params; i++)
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
  
  
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood(env, paramDomain, data_fname, &rxnMain);

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
  for (unsigned int i = 0; i < n_params; i++)
  {
      paramInitials[i] = init_kinetics[i];
  }

  // Now set up hypermeans of the Arrhenius inadequacy parameters
  for (unsigned int i = 0; i < n_params; i++)
  {
      paramInitials[n_params + i] = init_kinetics[i];
  }

  // Finally do the hypervariances
  for (unsigned int i = 0; i < n_params; i++)
  {
      paramInitials[2 * n_params + i] = 10.0;
  }

  // Set up covariance matrix for parameters
  QUESO::GslVector diagVec(paramSpace.zeroVector());
  QUESO::GslMatrix proposalCovMatrix(diagVec);

  double var = 1.0e-02;
  for (unsigned int i = 0; i < n_params_total; i++)
  {
      if (paramInitials[i] == 1.0e-16)
      {
         proposalCovMatrix(i,i) = var;
      }
      else
      {
         proposalCovMatrix(i,i) = var * paramInitials[i] * paramInitials[i];
      }
  }

  // Solve!
  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);
  
  return;
}
