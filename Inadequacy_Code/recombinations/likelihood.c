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
#include "model.h"
#include "truth_data.h"
// C functions 
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
//antioch
#include <antioch/kinetics_evaluator.h>
#include <antioch/kinetics_parsing.h>
// grvy
#include<grvy.h>
#include<sys/time.h>
#include<time.h>
// Eigen functions
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Set up necessary Eigen methods
using Eigen::EigenSolver;

template<class V, class M>
Likelihood<V, M>::Likelihood(const QUESO::BaseEnvironment& env, const QUESO::VectorSet<V, M> & domainSet, reaction_info * rxnInfo)
  : 
    QUESO::BaseScalarFunction<V, M>("", domainSet), 
    m_env(&env),                                                                   // QUESO Environment
    m_rxnMain(rxnInfo),                                                            // Problem and reaction information
    n_phis(rxnInfo->ProblemInfo.n_phis),                                           // Number of equivalence ratios
    n_times(rxnInfo->ProblemInfo.n_times),                                         // Number of sample points
    n_scen(rxnInfo->ProblemInfo.n_scenario),                                       // Number of scenarios (not including phi)
    num_fields(rxnInfo->ProblemInfo.n_species + rxnInfo->ProblemInfo.n_extra + 1), // Number of fields (species + temperature NOT including catchalls)
    n_species(rxnInfo->ProblemInfo.n_species + rxnInfo->ProblemInfo.n_atoms),      // Number of species (including catchalls)
    obs_data("truth_data.h5", n_phis, n_scen, n_times, num_fields)                 // Observation data class
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

  // Set up problem parameters
  const int n_atoms = rxn->ProblemInfo.n_atoms;             // Number of atoms in the model
  const int n_total = n_species + rxn->ProblemInfo.n_extra; // Species + catchalls + inert species (N2)
  const int n_data = rxn->ProblemInfo.n_species;            // Active species (not including inert species)
  const int dim = num_fields + n_atoms;                     // Total number of fields (catchalls, species, temperature)
  double varY   = 5.0e-03;
  double varT   = 2.0e+03;
  double var_ig = 10.0;

  // Set up initial conditions
  std::vector<double> initial_conditions(dim, 0.0);
  initial_conditions[1]  = rxn->ProblemInfo.oxidizer_i; // O2
  initial_conditions[7] = rxn->ProblemInfo.nitrogen;   // N2 

  // Solution vector
  std::vector<double> returnValues(n_times * dim + 1, 0.0); // Added one for ignition temperature

  // Set up inadequacy model parameters
  // First do Arrhenius parameters
  for (int j = 0; j < rxn->inad_model.n_reactions_inad; j++)
  {
      // Prefactor on j is 3 b/c of 3 Arrhenius params
      rxn->inad_model.Aj(j) = rxn->scale_factors[j] * exp(paramValues[3 * j]);
      rxn->inad_model.bj(j) = rxn->scale_factors[j] * paramValues[3 * j + 1];
      rxn->inad_model.Ej(j) = rxn->scale_factors[j] * paramValues[3 * j + 2];
  }

  // Then do thermo-chemistry parameters
  int n_arr_params = 3 * rxn->inad_model.n_reactions_inad;
  for (int k = 0; k < rxn->inad_model.n_atoms; k++)
  {
      for (int i = 0; i < 3; i++)
      {
          // i < 4 b/c using quadratics for thermo-chemistry and +1 for betas
          rxn->inad_model.alphas(k,i) = rxn->scale_factors[n_arr_params + 4 * k + i] * paramValues[n_arr_params + 4 * k + i];
      }
      rxn->inad_model.betas(k) = rxn->scale_factors[n_arr_params + 4 * k + 3] * paramValues[n_arr_params + 4 * k + 3];
  }

  double misfitValue = 0.0; // Difference between data and model
  double diff = 0.0;        // Argument of exponential in likelihood
  double fuel = rxn->ProblemInfo.fuel; // H2

  // Check if doing a heating rate problem.
  if (rxn->ProblemInfo.heat_rates)
  {
     initial_conditions[n_total] = rxn->ProblemInfo.TO;
  }

  std::vector<double> sample_points(n_times, 0.0);
  int scen; // For counting which scenario we're on

  grvy_timer_init("TIMING LIKELIHOOD EVALS"); // Initialize GRVY timer
  for (int i = 0; i < n_phis; ++i)
  { // Loop over equivalence ratio
      for (int ii = 0; ii < n_scen; ++ii)
      { // Loop over initial temperatures
          scen = i * n_scen + ii;
          initial_conditions[0] = fuel * obs_data.scenario_params[2*scen]; // Amount of H2
          // Check if doing a heating rate problem
          if (rxn->ProblemInfo.heat_rates)
          { // Set heating rate
             rxn->ProblemInfo.heating_rate = obs_data.scenario_params[2*scen+1]; 
          }
          else if (rxn->ProblemInfo.init_temperatures)
          { // Otherwise set initial temperature
             initial_conditions[n_total] = rxn->ProblemInfo.TO; // Temperature
          }
          // Get ignition data for this scenario
          rxn->ProblemInfo.time_ig = obs_data.ignition_data[2*scen];
          rxn->ProblemInfo.Tig     = obs_data.ignition_data[2*scen + 1];
          try
          { // Run the forward model
              for (int n = 0; n < n_times; n++)
              {
                  sample_points[n] = obs_data.sample_points[(i * n_scen + ii)* n_times + n];
              }
              if (scen > 0)
              {
                 grvy_timer_reset(); // Reset timer
              }
              grvy_timer_begin("Time forward model");
              hydrogenComputeModel(initial_conditions,sample_points,rxn,returnValues);
              grvy_timer_end("Time forward model");
              grvy_timer_finalize();
              grvy_timer_summarize();
              for (int j = 0; j < n_times; j++)
              { // Loop over the sample times
                  for (int k = 0; k < n_data; k++)
                  { // Loop over the species
                      // Calculate (d - model) where d is the data
                      diff = returnValues[dim * j + k] - 
                               obs_data.observation_data[n_times*num_fields*i + num_fields*j +k];
                      // Calculate (d - model)^T * \Sigma^{-1} * (d-model)
                      // Assume that \Sigma is a diagonal matrix where each entry is identical and equal to variance
                      misfitValue += diff * diff / varY;
                  }
                  // Temperature misfit
                  diff = returnValues[dim*j + dim - 1] - obs_data.observation_data[n_times*num_fields*i + num_fields*j + num_fields -1];
                  misfitValue += diff * diff / varT;
              }
              // Ignition temperature misfit
              diff = returnValues[n_times * dim] - rxn->ProblemInfo.Tig;
              misfitValue += diff * diff / var_ig;
          }
          catch(int exception)
          {
              misfitValue = 2.0e+30; // Get rid of solutions that caused and exception to be thrown.
          }
      } // End scenario loop
  } // End equivalence ratio loop

  std::cout << "   " << std::endl;
  std::cout << "log(L) = " << -0.5 * misfitValue << std::endl << std::endl;;

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

