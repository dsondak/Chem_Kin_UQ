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
  int n_phis;              // Equivalence ratios to run
  int n_heating;           // Different heating rates to run
  int n_T;                 // Different starting temperatures to run
  double timePoint;        // Time-step size
  int n_times;             // Number of time-steps to run
  int n_reactions;         // Number of reactions in mechanism
  double fuel;             // Stoichiometric factor for fuel (H2 here)
  double oxidizer_i;       // Initial concentration of oxidizer
  double nitrogen;         // Concentration of nitrogen
  int heat_rates_in;       // Flag for heating rate calibration
  int init_temperature_in; // Flag for initial temperature calibration
  double TO;               // Initial temperature
  double heating_rate;     // Heating rate
  char *thermo_filename;   // Thermodynamics input file
  char *reaction_filename; // Reaction input file
  double time_ig;          // Ignition time from detailed model
  double Tig;              // Ignition temperature from detailed model

  // Open input parameters file
  grvy_input_fopen("./input.txt");

  // Read in parameters
  grvy_input_fread_int("n_species", &n_species);
  grvy_input_fread_int("n_atoms", &n_atoms);
  grvy_input_fread_int("n_extra", &n_extra);
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
  grvy_input_fread_double("fuel", &fuel);
  grvy_input_fread_double("oxidizer_i", &oxidizer_i);
  grvy_input_fread_double("nitrogen", &nitrogen);
  grvy_input_fread_char("thermo", &thermo_filename);
  grvy_input_fread_char("reactionset", &reaction_filename);
  grvy_input_fread_double("time_ig", &time_ig);
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

  //scaling factors used for parameters
  // The last two are for the alpha parameters in the catchall reactions.
  const unsigned int n_scaling = n_ks + n_atoms;
  double scaling[n_scaling] = {
    1.0e10, 1.0e-2, 1.0e3, 1.0e7, 1.0e8,
    1.0e0,  1.0e0,  1.0e0, 1.0e0, 1.0e0, 
    1.0e3,  1.0e3,  1.0e3, 1.0e3, 1.0e3,
    1.0, 1.0};

  std::vector<double> scales(n_scaling,0.0);
  for (unsigned int i = 0; i < n_scaling; i++)
  {
     scales[i] = scaling[i];
  }

  // This class contains all the chemistry and problem size information
  reaction_info rxnMain(
      &chem_mixture,
      &reaction_set,
      &thermo,
      &kinetics,
      scales,
      n_species,
      n_atoms,
      n_extra,
      n_phis,
      n_scenario,
      n_times,
      n_reactions,
      oxidizer_i,
      nitrogen,
      fuel,
      heat_rates,
      init_temperatures,
      heating_rate,
      TO,
      time_ig,
      Tig);

  const unsigned int n_xi = rxnMain.S.nnz + rxnMain.S.n_eq_nonlin; // the linear operator and nonlinear catchall terms
  const unsigned int n_params = 3*n_xi + 9*n_atoms + 2;     // total number of parameters to be calibrated
                                                            // + 2 for the global Arrhenius parameter

  // Note:  The coefficient in front of n_xi is b/c we will calibrate xi but each
  //        xi has two hyperparameters to be calibrated.  The coefficient in front
  //        of n_atoms arises in the same manner only for the catchall reactions.

  std::cout << " " << std::endl;
  std::cout << "n_species = " << rxnMain.S.n_species << std::endl;
  std::cout << "n_atoms = " << rxnMain.S.n_atoms << std::endl;
  std::cout << "n_extra = " << rxnMain.S.n_extra << std::endl;
  std::cout << "n_eq_nonlin = " << rxnMain.S.n_eq_nonlin << std::endl;
  std::cout << "nnz = " << rxnMain.S.nnz << std::endl;
  std::cout << "n_ks = " << n_ks << std::endl;
  std::cout << "n_xi = " << n_xi << std::endl;
  std::cout << "n_params = " << n_params << std::endl;
  std::cout << "dim = " << dim << std::endl;
  std::cout << "n_phis = " << rxnMain.ProblemInfo.n_phis << std::endl;
  std::cout << "n_scenario = " << rxnMain.ProblemInfo.n_scenario << std::endl;
  std::cout << "n_times = " << rxnMain.ProblemInfo.n_times << std::endl;
  std::cout << "n_reactions = " << rxnMain.ProblemInfo.n_reactions << std::endl;
  std::cout << "oxidizer_i = " << rxnMain.ProblemInfo.oxidizer_i << std::endl;
  std::cout << "nitrogen = " << rxnMain.ProblemInfo.nitrogen << std::endl;
  std::cout << "fuel = " << rxnMain.ProblemInfo.fuel << std::endl;
  std::cout << "heating rate on = " << rxnMain.ProblemInfo.heat_rates << std::endl;
  std::cout << "init temperature on = " << rxnMain.ProblemInfo.init_temperatures << std::endl;
  std::cout << "heating rate = " << rxnMain.ProblemInfo.heating_rate << std::endl;
  std::cout << "initial temperature = " << rxnMain.ProblemInfo.TO << std::endl;
  std::cout << "ignition time (detailed)  = " << rxnMain.ProblemInfo.time_ig << std::endl;
  std::cout << "ignition temperature (detailed) = " << rxnMain.ProblemInfo.Tig << std::endl;
  std::cout << " " << std::endl;

  /*======================================================
  ***
  *** SIP Step 1 of 6: Instantiate the parameter space
  ***
  ========================================================*/
  
  QUESO::VectorSpace<QUESO::GslVector,QUESO::GslMatrix> paramSpace(env, "param_", n_params, NULL);
  
  /*======================================================
  ***
  *** SIP Step 2 of 6: Instantiate the parameter domain
  ***
  ========================================================*/

  QUESO::GslVector paramMinValues(paramSpace.zeroVector());
  QUESO::GslVector paramMaxValues(paramSpace.zeroVector());

  // The parameters have particular distributions.  The number of each parameter type
  // in paramInitials is different.  For example, \xi has a certain number of parameters
  // and that is different than the number of parameters in k_s.  So, the variables
  // below define the extent of the parameters in the array paramInitials.
  int mu_xi_range      = n_xi;
  int eta_xi_range     = mu_xi_range + n_xi;
  int xi_range         = eta_xi_range + n_xi;
  int mu_alpha0_range  = xi_range + n_atoms;
  int mu_alphai_range  = mu_alpha0_range + 2*n_atoms;
  int eta_alphai_range = mu_alphai_range + 3*n_atoms;
  int alpha0_range     = eta_alphai_range + n_atoms;
  int alphai_range     = alpha0_range + 2*n_atoms;
  int global_arr_range = alphai_range + 2;

  // Mean of xi
  for (int i = 0; i < mu_xi_range; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] = INFINITY;
  }
  
  // Variance of xi
  for (int i = mu_xi_range; i < eta_xi_range; i++)
  {
      paramMinValues[i] = 0.;
      paramMaxValues[i] = INFINITY;
  }

  // Xi
  for (int i = eta_xi_range; i < xi_range; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] = INFINITY;
  }

  // Mean of \alpha_{0k}
  for (int i = xi_range; i < mu_alpha0_range; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] = INFINITY;
  }

  // Mean of \alpha_{ik}, i > 0
  for (int i = mu_alpha0_range; i < mu_alphai_range; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] = INFINITY;
  }

  // Variance of \alpha_{ik}, i >= 0
  for (int i = mu_alphai_range; i < eta_alphai_range; i++)
  {
      paramMinValues[i] = 0.;
      paramMaxValues[i] = INFINITY;
  }

  // \alpha_{0k}
  for (int i = eta_alphai_range; i < alpha0_range; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] = INFINITY;
  }
  // \alpha_{ik}, i > 0
  for (int i = alpha0_range; i < alphai_range; i++)
  {
      paramMinValues[i] = -INFINITY;
      paramMaxValues[i] = INFINITY;
  }

  // Parameter range for global Arrhenius reaction
  paramMinValues[alphai_range] = -INFINITY;
  paramMaxValues[alphai_range] =  INFINITY;

  paramMinValues[alphai_range + 1] = -INFINITY;
  paramMaxValues[alphai_range + 1] =  INFINITY;

  QUESO::BoxSubset<QUESO::GslVector,QUESO::GslMatrix>
    paramDomain("param_", paramSpace, paramMinValues, paramMaxValues);

  /*======================================================
  ***
  *** SIP Step 3 of 6: Instantiate the likelihood function 
  *** object to be used by QUESO.
  ***
  ========================================================*/
  
  Likelihood<QUESO::GslVector, QUESO::GslMatrix> lhood(env, paramDomain, &rxnMain);

  /*======================================================
  ***
  *** SIP Step 4 of 6: Define the prior
  ***
  ========================================================*/

  // Chemistry prior for all parameters
  //QUESO::ChemistryVectorRV<QUESO::GslVector, QUESO::GslMatrix> priorTotal("priorTotal_", paramDomain);
  QUESO::ChemistryVectorRV<QUESO::GslVector, QUESO::GslMatrix> priorTotal("priorTotal_", paramDomain, n_xi, n_ks, n_atoms);

  /*======================================================
  ***
  *** SIP Step 5 of 6: Instantiate the inverse problem
  ***
  ========================================================*/

  QUESO::GenericVectorRV<>postTotal("post_", paramSpace); // Note extra prefix before the default "rv_" prefix
           
  QUESO::StatisticalInverseProblem<> ip("", NULL, priorTotal, lhood, postTotal); 
    
  /*======================================================
  ***
  *** SIP Step 6 of 6: Solve the inverse problem, that is,
  *** set the 'pdf' and the 'realizer' of the posterior RV
  ***
  ========================================================*/
  std::cout << "Solving the SIP with Multi-Level Metropolis Hastings" << std::endl << std::endl;  
 
  //The following is set if using ip.solveWithBayesMetropolisHastings
  QUESO::GslVector paramInitials(paramSpace.zeroVector());

  // Initialize parameter values
  double init_xi = 3.0;
  double init_alpha0 = 10.0;
  double init_alphai = -6.0;
  // Mean of the prior mean for xi and the catchall reactions
  for (int i = 0; i < mu_xi_range; i++)
  {
      paramInitials[i] = init_xi;
  }
  // Variance of the prior variance for xi and the catchall reactions
  for (int i = mu_xi_range; i < eta_xi_range; i++)
  {
      paramInitials[i] = 80.0;
  }
  // xi and the catchall reactions
  for (int i = eta_xi_range; i < xi_range; i++)
  {
      paramInitials[i] = init_xi;
  }
  // Catchall enthalpies:
  //   hk = \alpha_{0k} + \alpha_{1k}*T + \alpha_{2k} * T * T
  //
  // Mean of the prior mean for \alpha_{0k}
  for (int i = xi_range; i < mu_alpha0_range; i++)
  {
      paramInitials[i] = init_alpha0;
  }
  // Mean of the prior mean for \alpha_{ik}, i > 0
  for (int i = mu_alpha0_range; i < mu_alphai_range; i++)
  {
      paramInitials[i] = init_alphai;
  }
  // Variance of the prior variance for all catchall enthalpy coefficients
  for (int i = mu_alphai_range; i < eta_alphai_range; i++)
  {
      paramInitials[i] = 80.0;
  }
  // Catchall enthalpy coefficients for \alpha_{0k}
  for (int i = eta_alphai_range; i < alpha0_range; i++)
  {
      paramInitials[i] = init_alpha0;
  }
  // Catchall enthalpy coefficients for \alpha_{ik}, i > 0
  for (int i = alpha0_range; i < alphai_range; i++)
  {
      paramInitials[i] = init_alphai;
  }
  // Initial value of global Arrhenius activation energy
  //paramInitials[alphai_range    ] = 1000.0;
  paramInitials[alphai_range] = 919.0;
  paramInitials[alphai_range + 1] = 2680.0;

  // Set up covariance matrix for parameters
  QUESO::GslVector diagVec(paramSpace.zeroVector());
  QUESO::GslMatrix proposalCovMatrix(diagVec);

  double var = 1.0e-04;
  for (int i = 0; i < n_params; i++)
  {
      proposalCovMatrix(i,i) = var;
  }
  proposalCovMatrix(n_params-2, n_params-2) = 8445.61;
  proposalCovMatrix(n_params-1, n_params-1) = 71824.0;

  ip.solveWithBayesMetropolisHastings(NULL, paramInitials, &proposalCovMatrix);

  return;
}
