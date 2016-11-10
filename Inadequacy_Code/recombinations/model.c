/*-------------------------------------------------------------------
 *
 ***   model.c   ***
 *
 * This is the code for the forward model of the 0D reactor.
 *
 * This file contains three functions:
 *     1.)  hydrogenComputeModel
 *     2.)  hydrogenFunction
 *     3.)  hydrogenJacobian
 *
 * 1.) hydrogenComputeModel
 *     INPUTS:
 *            o initialValues: Initial conditions for each species
 *                             and temperature
 *            o timePoints   : Times at which to take data
 *            o rxn         : Class containing all information about
 *                             chemical reactions, thermodynamics, and
 *                             problem size.
 *
 *     OUTPUTS:
 *             o returnValues: Species and temperature at each of the
 *                             times to be sampled.
 *
 * 2.) hydrogenFunction
 *     INPUTS:
 *            o t: The time
 *            o Y: Vector of solutions (species and temperature)
 *            o dYdt: Right-hand-sides of each ODE
 *            o params: Parameters including chemistry, thermodynamics
 *                      and problem size
 *
 *     OUTPUTS:
 *             o Y: Solutions at each time of interest
 *
 * Contact:  David Sondak, dsondak@ices.utexas.edu
 *
 *-----------------------------------------------------------------*/
// User-defined functions
#include "model.h"
#include "reaction_info.h"
// C functions
#include <cmath>
#include <vector>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <algorithm>
// GSL functions
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
// Eigen functions
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
//antioch
#include <antioch/kinetics_evaluator.h>
#include <antioch/cea_evaluator.h>

#ifndef __EPS_ABS
#define __EPS_ABS 1e-8
#endif
#ifndef __EPS_REL
#define __EPS_REL 1e-8
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;

/***************************
 *
 * hydrogenFunction
 * 
 ***************************/
int hydrogenFunction(double t, const double Y[], double dYdt[], void* params)
{
  // Unpack parameters
  reaction_info rxn = *(reaction_info *) params;

  // Get the number of species in the model
  unsigned int n_species_total = rxn.ProblemInfo.n_species + rxn.ProblemInfo.n_extra;
  unsigned int n_atoms   = rxn.ProblemInfo.n_atoms;
  int dim = n_species_total + n_atoms + 1;

  // Solution vector (units are moles)
  std::vector<double> molar_densities(n_species_total, 0.0);
  for (unsigned int i = 0; i < n_species_total; i++)
  { // Populate molar_densities vector
      molar_densities[i] = Y[i];
      if(molar_densities[i] <= 0)
      { // Guard against negative concentrations
        molar_densities[i] = 0;
      }
  }

  // Prepare Antioch to get chemistry and thermodynamic information
  double temperature = Y[dim-1];
  Antioch::TempCache<double> temp_cache(temperature);

  // Get some thermodynamic data for each species
  std::vector<double> h_RT_minus_s_R(n_species_total);
  rxn.Thermo->h_RT_minus_s_R(temp_cache, h_RT_minus_s_R);

  // Perform kinetics calculations 
  // (basically just gives RHS of ODEs for each species)
  // (Note that temperature is passed in but the RHS for the
  //  temperature equation is NOT computed by Antioch.  Just
  //  need temperature in reaction rate calculations.)
  std::vector<double> mole_sources(n_species_total, 0.0);
  rxn.Kinetics->compute_mole_sources(
      temperature,     // temperature needed for calculations
      molar_densities, // current concentrations (in moles)
      h_RT_minus_s_R,  // exponent in equilibrium constant for reverse rates
      mole_sources);   // RHS

  // Set up RHS of each ODE for the species
  for (unsigned int i = 0; i < n_species_total; i++)
  { // Loop over species
      dYdt[i] = mole_sources[i];
  }

  // Needed for energy equation
  double Qnum = 0.0;
  double Qden = 0.0;

  double R_universal = Antioch::Constants::R_universal<double>();
  double RT = R_universal * temperature;

  // Get number of species in inad. model.  Here it's just 
  // the usual species minus N2, H and O plus H' and O'
  const unsigned int n_species_inad = n_species_total  - 3 + n_atoms;

  // Create a new vector containing only relevant species
  std::vector<double> Yinad(n_species_inad, 0.0);

  // Skip H and O
  for (int k = 0; k < 2; k++)
  {
      Yinad[k] = Y[k];
  }
  // Skipped H and O so decrement index by 2
  for (int k = 4; k < rxn.ProblemInfo.n_species; k++)
  {
      Yinad[k - 2] = Y[k];
  }
  // Next copy catchall species to Yinad
  for (int k = 0; k < n_atoms; k++)
  {
      Yinad[rxn.ProblemInfo.n_species + k - 2] = Y[(dim - 1) - n_atoms + k];
  }

  // Compute catchall enthalpy and specific heat
  rxn.inad_model.calc_h_prime(temperature);
  rxn.inad_model.calc_cp_prime(temperature);

  // Before computing catchall entropy, we need to 
  // compute the entropy coefficients based off of the 
  // constraints.  Let's do that right now.

  double g1 = h_RT_minus_s_R[0]; // H2
  double g2 = h_RT_minus_s_R[1]; // O2
  double g3 = h_RT_minus_s_R[4]; // OH
  double g4 = h_RT_minus_s_R[5]; // HO2
  double g5 = h_RT_minus_s_R[6]; // H2O

  double s1T;
  double s2T;

  s1T = rxn.inad_model.alphas(0,1) * log(temperature) + 2.0 * rxn.inad_model.alphas(0,2) * temperature;
  s2T = rxn.inad_model.alphas(1,1) * log(temperature) + 2.0 * rxn.inad_model.alphas(1,2) * temperature;

  std::vector<double> constraints(5, 0.0); // Vector of constraints

  // First constraint
  constraints[0] = -s1T + (rxn.inad_model.h_prime(0) - 0.5 * g1) / temperature;

  // Directly impose constraint on \beta_1
  rxn.inad_model.betas(0) = constraints[0] - rxn.inad_model.betas(0);

  // Now we need to calculate the remaining constraints and determine which one 
  // is the most strict
  constraints[0] = -s2T + (rxn.inad_model.h_prime(1) - 0.5 * g2) / temperature;
  constraints[1] = -rxn.inad_model.betas(0) - s1T - s2T + 
          (rxn.inad_model.h_prime(0) + rxn.inad_model.h_prime(1) - g3) / temperature;
  constraints[2] = 0.5 * (-rxn.inad_model.betas(0) - s1T - 2.0 * s2T + 
                (rxn.inad_model.h_prime(0) + 2.0 * rxn.inad_model.h_prime(1) - g4) / temperature);
  constraints[3] = -2.0 * rxn.inad_model.betas(0) - 2.0 * s1T - s2T + 
          (2.0 * rxn.inad_model.h_prime(0) + rxn.inad_model.h_prime(1) - g5) / temperature;

  double min_constraint = *std::min_element(constraints.begin(), constraints.end());

  // Now we're ready to assign the second entropy coefficient
  rxn.inad_model.betas(1) = min_constraint - rxn.inad_model.betas(1);

  // And now compute the entropy
  rxn.inad_model.calc_s_prime(temperature);

  // Set up s/R - h_RT for each species
  std::vector<double> delta_k(n_species_inad, 0.0);

  // First just do H2 and O2
  for (int k = 0; k < 2; k++)
  {   
      delta_k[k] = rxn.Thermo->s_over_R(temp_cache, k) - 
                   rxn.Thermo->h_over_RT(temp_cache, k);
  }
  // Skipped H and O so need to decrement index by 2
  for (int k = 4; k < rxn.ProblemInfo.n_species; k++)
  {   
      delta_k[k - 2] = rxn.Thermo->s_over_R(temp_cache, k) - 
                   rxn.Thermo->h_over_RT(temp_cache, k);
  }
  // Now do catchalls (still need to decrement by 2)
  for (int k = 0; k < n_atoms; k++)
  {
      delta_k[rxn.ProblemInfo.n_species + k - 2] = rxn.inad_model.s_prime(k) / R_universal - 
                                                   rxn.inad_model.h_prime(k) / RT;
  }

  rxn.inad_model.progress_rate(Yinad, temperature, R_universal, delta_k);

  // Finally, compute the RHS for these species.
  VectorXd omega_dot_inad;
  omega_dot_inad = rxn.inad_model.nukj * rxn.inad_model.rj;

  // The very last step is to add omega_dot_inad to the 
  // dYdt and then move on to the energy equation.
  // First just do H2 and O2
  for (int k = 0; k < 2; k++)
  {
      dYdt[k] += omega_dot_inad(k);
  }
  // Skipping H and O so decrement index by 2
  for (int k = 4; k < rxn.ProblemInfo.n_species; k++)
  {
      dYdt[k] += omega_dot_inad(k - 2);
  }
  // Now do catchalls
  for (int k = 0; k < n_atoms; k++)
  {
      dYdt[(dim - 1) - n_atoms + k] = omega_dot_inad(rxn.ProblemInfo.n_species + k - 2);
  }

  // Get catchall contribution to energy equation
  for (unsigned int k = 0; k < n_atoms; k++)
  {
      Qnum += rxn.inad_model.h_prime(k) * dYdt[(dim - 1) - n_atoms + k];
      Qden += rxn.inad_model.cp_prime(k) * Y[(dim - 1) - n_atoms + k];
  }

  // Energy equation computations
  double enthalpy; // Enthalpy
  double specific_heat_p;  // cp

  for (unsigned int s = 0; s < n_species_total; s++)
  { // Get numerator and denominator in energy equation (sum of species)
      // Get enthalpy and convert to molar from mass basis
      // Note that R_universal/Rs = Ws where Ws is the molecular weight of species s
      enthalpy = rxn.Thermo->h(temp_cache,s) * R_universal / rxn.Chem_mixture->R(s);
      // get cp and convert from mass to molar
      specific_heat_p = rxn.Thermo->cp(temp_cache,s) * R_universal / rxn.Chem_mixture->R(s);
      // Numerator in energy equation: h_s * dx_s/dt
      Qnum  += enthalpy * mole_sources[s];
      // Denominator in energy equation: cp_s * x_s
      Qden += specific_heat_p * molar_densities[s];
  }

  if (Qden <= 0)
  {
      std::cout << "Disaster!  Qden <=0.  In fact, Qden = " << Qden << std::endl;
      exit(0);
  }

  // Right hand side of energy equation.
  dYdt[dim-1] = (-Qnum + rxn.ProblemInfo.heating_rate) / Qden;
  
  return GSL_SUCCESS;
}

/***************************
 *
 * hydrogenJacobian
 * 
 ***************************/
int hydrogenJacobian(double t, const double Y[], double *dfdY, double dfdt[], void* params)
{
  return GSL_SUCCESS;
}

/***************************
 *
 * hydrogenComputeModel
 * 
 ***************************/
void hydrogenComputeModel(
  std::vector<double>&  initialValues,
  std::vector<double>&  timePoints,
  reaction_info*        rxn,  
  std::vector<double>&  returnValues)
{  
  unsigned int dim = initialValues.size(); // Size of system

  // Set up GSL solver
  gsl_odeiv2_system sys = {hydrogenFunction, hydrogenJacobian, dim, rxn};
  
  // Pass parameters to GSL ODE solver
  double dt      = 1.0e-06; // initial time step size
  double err_abs = 1.0e-08;  // Absolute error tolerance for adaptive stepping
  double err_rel = 1.0e-08;   // Relative error tolerance for adaptive stepping
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, dt, err_abs, err_rel);   

  // Initialize solution field
  double Y[dim];
  for (unsigned int i = 0; i < dim; ++i)
  {
      Y[i] = initialValues[i];
  }

  // Prepare for time integrator to get ignition data
  double t = 0.0; // Initial time
  double dt_samp = 1e-6;
  int Nt = 100000;
  double finalTime;

  bool ignition = false;

  double Tig = 0.0;
  double time_ig = 0.0;

  int status;

  // Time integration
  for (int i = 0; i < Nt; i++)
  { // Integration loop
      finalTime = i * dt_samp; // Integration time
      while (t < finalTime)
      {
         // Call GSL ODE integrator 
         // t is the current time
         // finalTime is the integration time
         // Y is the solution at time t
         status = gsl_odeiv2_driver_apply( d, &t, finalTime, Y );
         if (status != GSL_SUCCESS)
         {
            std::cout << "GSL wants to a take dt < dt_min." << std::endl;
            throw status;
         }
         if (Y[0] < 0.999*initialValues[0])
         {
            Tig = Y[dim-1];
            time_ig = t;
            ignition = true;
         }
      }

      if (ignition)
      {
         break;
      }

  } // end loop over time points

/*=========================================
 *
 * Determine time shift
 *
 *=========================================*/
  // New version:  Makes sure the shifted time
  // does not become negative.  Shifts time
  // so that the ignition times coincide.
  // Specifies the second point arbitrarily
  // to enforce non-negativity.

  double scale;
  double shift;
  double delta_tig = rxn->ProblemInfo.time_ig - time_ig;

  int n_samp = rxn->ProblemInfo.n_times;
  double t_fit;
  t_fit = - timePoints[n_samp - 1] * time_ig / (delta_tig + timePoints[n_samp - 1]);
  scale = (time_ig - t_fit) / (rxn->ProblemInfo.time_ig - t_fit);
  shift = t_fit * delta_tig / (rxn->ProblemInfo.time_ig - t_fit);

  std::vector<double> timePoints_shift(n_samp, 0.0);
  for (int i = 0; i < n_samp; i++)
  {
      timePoints_shift[i] = scale * timePoints[i] + shift;
      if (timePoints_shift[i] < 0)
      {
         std::cout << "The time shift failed" << std::endl;
         std::cout << "b/c it tried to make " << std::endl;
         std::cout << "time negative.       " << std::endl;
         exit(0);
      }
  }

  // Reinitialize solution field
  for (int i = 0; i < dim; i++)
  {
      Y[i] = initialValues[i];
  }

  // Prepare for time integrator
  t = 0.0; // Initial time

  // Time integration
  for (unsigned int i = 0; i < n_samp; i++)
  { // Integration loop
      //finalTime = timePoints_shift[i];
      finalTime = timePoints[i];
      // Call GSL ODE integrator 
      // t is the current time
      // finalTime is the integration time
      // Y is the solution at time t
      status = gsl_odeiv2_driver_apply( d, &t, finalTime, Y );
      if (status != GSL_SUCCESS)
      {
         std::cout << "GSL wants to take dt < dt_min." << std::endl;
         throw status;
      }
      //std::cout << "Time = " << t << std::endl;

      // Store results in return field
      for (unsigned int j = 0; j < dim; j++)
      {
          returnValues[dim * i + j] = Y[j];
      }

  } // end loop over time points
  // Don't forget to store the ignition temperature!
  returnValues[dim * n_samp] = Tig;

  // deallocate memory   
  gsl_odeiv2_driver_free( d );

} // end hydrogenComputeModel
