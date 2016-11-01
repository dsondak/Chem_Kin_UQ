/*-------------------------------------------------------------------
 *
 ***   qoi.c   ***
 *
 * This script defines the qoi.
 *
 *     INPUTS:
 *             o env:       QUESO environment
 *             o phis:      Vector of equivalence ratios
 *             o scenario:  Vector of initial temperatures
 *                          or heating rates
 *             o times:     Sample times
 *             o concs:     Species concentrations data
 *             o rxnInfo:   Structure containing chemistry
 *                          and problem information
 *
 *     OUTPUTS:
 *             o qoiValues:  Computed QOI
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/

// User-defined functions
#include <qoi.h>
#include <model.h>
#include "reaction_info.h"
// C functions 
#include <cmath>
//antioch
#include <antioch/kinetics_evaluator.h>
#include <antioch/kinetics_parsing.h>

// Constructor for the qoi routine
qoiRoutine_Data::qoiRoutine_Data(const QUESO::BaseEnvironment& env,
    const std::vector<double> & phis,
    const std::vector<double> & scenario,
    const std::vector<double> & times,
    std::vector<double> & concs,
    reaction_info * rxnInfo)
: m_env(&env),
  m_phis(phis),
  m_scenario(scenario),
  m_times(times),
  m_concs(concs),
  m_rxnMain(rxnInfo)
{
}

// Destructor
qoiRoutine_Data::~qoiRoutine_Data()
{
}

/*========================================
 * 
 * The user-defined qoi routine
 * 
 ========================================*/

void qoiRoutine(
  const QUESO::GslVector&                    paramValues,
  const QUESO::GslVector*                    paramDirection,
  const void*                                functionDataPtr,
        QUESO::GslVector&                    qoiValues,
        QUESO::DistArray<QUESO::GslVector*>* gradVectors,
        QUESO::DistArray<QUESO::GslMatrix*>* hessianMatrices,
        QUESO::DistArray<QUESO::GslVector*>* hessianEffects)
{
  const QUESO::BaseEnvironment& env = paramValues.env();

  if (paramDirection && gradVectors && hessianEffects && hessianMatrices)
  {
    // Logic just to avoid warnings from INTEL compiler
  }
  
  // Unpack data from qoi data structure
  const std::vector<double>& phis = ((qoiRoutine_Data *) functionDataPtr)->m_phis; // Equivalence ratios
  const std::vector<double>& scenario = ((qoiRoutine_Data *) functionDataPtr)->m_scenario; // Initial temeratures
  const std::vector<double>& times = ((qoiRoutine_Data *) functionDataPtr)->m_times; // Sampel times
  const std::vector<double>& concs = ((qoiRoutine_Data *) functionDataPtr)->m_concs; // Species concentrations
  reaction_info *            rxn   = ((qoiRoutine_Data *) functionDataPtr)->m_rxnMain; // Reaction information and problem size


  // Set up problem parameters
  const int n_atoms = rxn->S.n_atoms;                   // Number of atoms in the model
  const int n_species = rxn->S.n_species + n_atoms;     // Number of species in the model (+n_atoms to account for catchalls)
  const int n_total = n_species + rxn->S.n_extra;       // Species + catchalls + neglected species (H2O2 and N2)
  const int dim = n_total + 1;                          // All fields (+1 for temperature)
  const int n_phis = rxn->ProblemInfo.n_phis;           // Number of equivalence ratios
  const int n_scenario = rxn->ProblemInfo.n_scenario;             // Number of temperatures
  const int n_times = rxn->ProblemInfo.n_times;         // Number of sample points
  const int n_reactions = rxn->ProblemInfo.n_reactions; // Number of reactions in the reduced model
  const int n_ks = 3 * n_reactions;                     // Number of parameters in the reduced model
  const int n_xi = rxn->S.nnz + rxn->S.n_eq_nonlin;     // Numer of parameters in the stochastic operator

  // Allocate equivalence ratio data
  std::vector<double> phiPoints(n_phis,0.);
  for (int j = 0; j < n_phis; j++)
  {
      phiPoints[j] = phis[n_times * n_scenario * j];
  }

  // Allocate temperature data
  std::vector<double> scenarioPoints(n_scenario,0.);
  for (int j = 0; j < n_scenario; j++)
  {
      scenarioPoints[j] = scenario[n_times * j];
  }

  // Allocate sample points
  std::vector<double> timePoints(n_times,0.);
  for (int i = 0; i < n_times; i++)
  {
      timePoints[i] = times[i];
  }

  // Set up initial conditions
  std::vector<double> initial_conditions(dim, 0.0);
  initial_conditions[3]  = rxn->ProblemInfo.oxidizer_i; //O2
  initial_conditions[10] = rxn->ProblemInfo.nitrogen;   //N2 

  // Solution vector
  std::vector<double> returnValues(n_times * dim,0.0);

  // Multiply newValues for A by constants to keep params close to 1
  std::vector<double> tempValues(n_ks);
  //GAUSSIAN
   // Transform gaussian to lognormal and scale
  for (int i = 0; i < n_reactions; i++) tempValues[i] = rxn->Scales[i] * exp(paramValues[i + 3*n_xi]);
  // Multiply newValues for A by constants to keep params close to 1
  for (int i = n_reactions; i < n_ks; i++) tempValues[i] = rxn->Scales[i] * paramValues[i + 3*n_xi];
  //END GAUSSIAN

  // Copy paramValues into struct
  rxn->model_params.resize(n_xi + rxn->S.n_eq_nonlin * rxn->S.n_atoms);
  for (int i = 0; i < n_xi; i++)
  {
      rxn->model_params[i] = paramValues[i + 2 * n_xi];
  }

  // Copy energies from paramValues into struct
  for (int i = 0; i < n_atoms; i++) 
    rxn->model_params[n_xi + i] = paramValues[3*n_xi + n_ks + 6*n_atoms + i] * rxn->Scales[n_ks + i];
  for (int i = n_atoms; i < 3*n_atoms; i++) 
    rxn->model_params[n_xi + i] = paramValues[3*n_xi + n_ks + 6*n_atoms + i];

  // Reset reaction rate parameters to rescale values
  std::vector<double> newValues(3);
  for (int i = 0; i < n_reactions; i++)
  {
      newValues[0] = tempValues[i];
      newValues[1] = tempValues[i + n_reactions];
      newValues[2] = tempValues[i + 2*n_reactions];
      reset_rate(rxn->Reaction_set->reaction(i).forward_rate(),newValues);
  }

  double fuel = rxn->ProblemInfo.fuel; // H2
  // Check if doing a heating rate problem

  if (rxn->ProblemInfo.heat_rates)
  {
     initial_conditions[n_total] = rxn->ProblemInfo.TO;
  }

  for (int i = 0; i < n_phis; i++)
  { // Loop over equivalence ratio
      for (int ii = 0; ii < n_scenario; ii++)
      { // Loop over initial temperatures
          initial_conditions[n_atoms] = fuel * phiPoints[i]; //amount of H2
          // Check if doing a heating rate problem
          if (rxn->ProblemInfo.heat_rates)
          { // Set heating rate
             rxn->ProblemInfo.heating_rate = scenarioPoints[ii];
          }
          else if (rxn->ProblemInfo.init_temperatures)
          { // Otherwise set initial temperature
             initial_conditions[n_total] = scenarioPoints[ii];        // Temperature
          }
          // Run the forward model
          hydrogenComputeModel(initial_conditions,timePoints,rxn,returnValues);
          for (int j = 0; j < returnValues.size(); j++)
          { // Return solution
              qoiValues[(n_scenario*n_times*dim) * i + (n_times*dim) * ii + j] = returnValues[j];
          }
      }
  }
  return;
}
