/*-------------------------------------------------------------------
 *
 ***   likelihood.c   ***
 *
 * This script computes the likelihood.
 *
 *     INPUTS:
 *             o env:       QUESO environment
 *             o phis:      Vector of equivalence ratios
 *             o scenario:  Vector of initial temperatures
 *                          or heating rates
 *             o times:     Sample times
 *             o concs:     Species concentrations data
 *             o temps:     Temperature data
 *             o rxnInfo:   Structure containing chemistry
 *                          and problem information
 *
 *     OUTPUTS:
 *             o misfitValue:  Logarithm of likelihood
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

// User-defined functions
#include "likelihood.h"
#include "reaction_info.h"
#include "inadequacy_model.h"
#include "kinetics_solve.h"
#include "truth_data.h"
// C functions 
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
//antioch
#include <antioch/kinetics_evaluator.h>
#include <antioch/kinetics_parsing.h>
#include <antioch/nasa_mixture.h>
#include <antioch/nasa_mixture_parsing.h>
#include <antioch/nasa_curve_fit_base.h>
// Eigen functions
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Set up necessary Eigen methods
using Eigen::EigenSolver;

template<class V, class M>
Likelihood<V, M>::Likelihood(const QUESO::BaseEnvironment& env, const QUESO::VectorSet<V, M> & domainSet, const char* fname, reaction_info * rxnInfo)
//Likelihood<V, M>::Likelihood(const QUESO::BaseEnvironment& env, const QUESO::VectorSet<V, M> & domainSet, reaction_info * rxnInfo)
  : 
    QUESO::BaseScalarFunction<V, M>("", domainSet), 
    m_env(&env),                                                                   // QUESO Environment
    m_rxnMain(rxnInfo),                                                            // Problem and reaction information
    n_phis(rxnInfo->ProblemInfo.n_phis),                                           // Number of equivalence ratios
    n_times(rxnInfo->ProblemInfo.n_times),                                         // Number of sample points
    n_scen(rxnInfo->ProblemInfo.n_scenario),                                       // Number of scenarios (not including phi)
    num_fields(rxnInfo->ProblemInfo.n_species_d + rxnInfo->ProblemInfo.n_inert + 1), // Number of fields (species + temperature NOT including catchalls)
    n_species(rxnInfo->ProblemInfo.n_species),                                     // Number of species (including catchalls)
    //obs_data("truth_data.h5", n_phis, n_scen, n_times, num_fields)                 // Observation data class
    obs_data(fname, n_phis, n_scen, n_times, num_fields)                 // Observation data class
{
  // Constructor
}

template <class V, class M>
Likelihood<V, M>::~Likelihood()
{
  // Deconstruct here
}

/*========================================
 * 
 * The user-defined likelihood routine
 * 
 ========================================*/

template<class V, class M>
double Likelihood<V, M>::lnValue(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect) const
{
  if (paramDirection && gradVector && hessianMatrix && hessianEffect) 
  {
    // Just to eliminate INTEL compiler warnings
  }
  
  // Unpack data from likelihood data structure
  reaction_info * rxn = m_rxnMain; // Reaction information and problem size

  //std::ofstream param_file = write_file;

  // Set up problem parameters
  const unsigned int n_inad  = rxn->InadInfo.n_inad;               // Number of atoms in the model
  const unsigned int n_data  = rxn->ProblemInfo.n_species;         // Active species (not including inert species)
  const unsigned int n_eq    = rxn->ProblemInfo.n_eq;              // Total number of fields (catchalls, species, temperature)
  const unsigned int n_inert = rxn->ProblemInfo.n_inert;           // Number of inert species (e.g. N2)
  double varY   = 5.0e-03;
  double varT   = 2.0e+03;
  double var_ig = 10.0;

  // Set up initial conditions
  std::vector<double> initial_conditions(n_eq, 0.0);
  initial_conditions[1] = rxn->ProblemInfo.oxidizer_i; // O2
  initial_conditions[n_eq-1] = rxn->ProblemInfo.nitrogen;   // N2 

  // Solution vector
  std::vector<double> returnValues(n_times * n_eq + 1, 0.0); // Added one for ignition time

  /* Set up inadequacy model parameters */
  // First reset Arrhenius parameters
  std::vector<double> newValues(3, 0.0);
  unsigned int n_reactions = rxn->ProblemInfo.n_reactions;
  unsigned int n_reactions_inad = rxn->InadInfo.n_reactions_inad;
  printf("Parameter Values:\n");
  for (int j = 0; j < n_reactions_inad; j++)
  {
      // Prefactor on j is 3 b/c of 3 Arrhenius params
      newValues[0] = rxn->scale_factors[3 * j    ] * exp(paramValues[3 * j]);
      newValues[1] = rxn->scale_factors[3 * j + 1] * paramValues[3 * j + 1];
      newValues[2] = rxn->scale_factors[3 * j + 2] * paramValues[3 * j + 2];
      reset_rate(rxn->Reaction_set->reaction(n_reactions + j).forward_rate(),newValues);
      printf("A = %25.16e     b = %25.16e     E = %25.16e\n", newValues[0], newValues[1], newValues[2]);
  }
  printf("\n\n");

  double misfitValue = 0.0;                   // Difference between data and model
  double diff        = 0.0;                   // Argument of exponential in likelihood
  double fuel        = rxn->ProblemInfo.fuel; // H2

  std::vector<double> sample_points(n_times, 0.0);
  int scen; // For counting which scenario we're on

  for (int i = 0; i < n_phis; ++i)
  { // Loop over equivalence ratio
      for (int ii = 0; ii < n_scen; ++ii)
      { // Loop over initial temperatures
          misfitValue = 0.0; // Difference between data and model
          scen = i * n_scen + ii;
          initial_conditions[0] = fuel * obs_data.scenario_params[2*scen]; // Amount of H2
          // Check if doing a heating rate problem
          if (rxn->ProblemInfo.heat_rates)
          { // Set heating rate
             rxn->ProblemInfo.heating_rate = obs_data.scenario_params[2*scen+1]; 
          }
          else if (rxn->ProblemInfo.init_temperatures)
          { // Otherwise set initial temperature
             initial_conditions[n_eq-1] = rxn->ProblemInfo.TO; // Temperature
          }
          // Get ignition data for this scenario
          rxn->ProblemInfo.time_ig = obs_data.ignition_data[2*scen];
          rxn->ProblemInfo.Tig     = obs_data.ignition_data[2*scen + 1];
          for (int n = 0; n < n_times; n++)
          {
              sample_points[n] = obs_data.sample_points[(i * n_scen + ii)* n_times + n];
          }
          try
          { // Run the forward model
              // Compute the solution
              kinetics_forward(initial_conditions,sample_points,rxn,returnValues);
              // Compute the misfit
              for (int j = 0; j < n_times; j++)
              { // Loop over the sample times
                  for (int k = 0; k < n_data; k++)
                  { // Loop over the species
                      // Calculate (d - model) where d is the data
                      diff = returnValues[n_eq * j + k] - 
                               obs_data.observation_data[n_times*num_fields*i + num_fields*j +k];
                      // Calculate (d - model)^T * \Sigma^{-1} * (d-model)
                      // Assume that \Sigma is a diagonal matrix where each entry is identical and equal to variance
                      misfitValue += diff * diff / varY;
                  }
              }
              // Ignition temperature misfit
              diff = returnValues[n_times * n_eq] - rxn->ProblemInfo.Tig;
              misfitValue += diff * diff / var_ig;
          }
          catch(...)
          {
              std::cout << "Faulty Parameters:\n\n";
              for (unsigned int j = 0; j < n_reactions_inad; j++)
              {
                  // Prefactor on j is 3 b/c of 3 Arrhenius params
                  std::cout << rxn->scale_factors[3 * j    ] * exp(paramValues[3 * j]) << "   "
                            << rxn->scale_factors[3 * j + 1] * paramValues[3 * j + 1]  << "   "
                            << rxn->scale_factors[3 * j + 2] * paramValues[3 * j + 2] <<"\n";
              }
              
              std::cout << "\n\n" << std::endl;
              misfitValue = 2.0e+30; // Get rid of solutions that caused and exception to be thrown.
          }
      } // End scenario loop
  } // End equivalence ratio loop

  printf("\n");
  printf("log(L) = %25.16e\n\n", -0.5 * misfitValue);

  return (-0.5 * misfitValue);  // Return likelihood

} // End lnValue


template<class V, class M>
double Likelihood<V, M>::actualValue(
  const QUESO::GslVector& paramValues,
  const QUESO::GslVector* paramDirection,
  QUESO::GslVector*       gradVector,
  QUESO::GslMatrix*       hessianMatrix,
  QUESO::GslVector*       hessianEffect) const
{
  return 0.0;
}

template class Likelihood<QUESO::GslVector, QUESO::GslMatrix>;

