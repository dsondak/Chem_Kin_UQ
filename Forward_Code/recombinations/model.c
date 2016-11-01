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

#include "model.h"
#include "reaction_info.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <assert.h>
#include <Eigen/Dense>

#ifndef __EPS_ABS
#define __EPS_ABS 1e-8
#endif
#ifndef __EPS_REL
#define __EPS_REL 1e-8
#endif

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::MatrixXi;


Eigen::VectorXd reaction_rate(std::vector<double> Yinad, double temperature, 
                              std::vector<double> h_prime, std::vector<double> s_prime, 
                              reaction_info rxn)
{
     double R_universal = Antioch::Constants::R_universal<double>();
     Antioch::TempCache<double> temp_cache(temperature);

     double RT = R_universal * temperature;
     double pa = 1.0e+05; // 1 bar in Pa
     double pa_RT = pa / RT;

     int n_atoms = rxn.ProblemInfo.n_atoms;

     MatrixXd nukj_r = rxn.ProblemInfo.nukj_r;
     MatrixXd nukj_p = rxn.ProblemInfo.nukj_p;
     MatrixXd nukj   = rxn.ProblemInfo.nukj;
     VectorXd gamma = rxn.ProblemInfo.gamma;

     int Mc = 5; // FIXME  Remove hard coding
     int Nc = rxn.ProblemInfo.n_species;
     double exp_arg_j;

     // Make up some Arrhenius coeffs. for now
     VectorXd A(Mc);
     VectorXd beta(Mc);
     VectorXd Ea(Mc);

     // 5 RXN TESTS 
     //A    << 1.0e+08, 2.0e-02, 1.0e+07, 1.5e+03, 1.17e+03;
     beta << 0.5, 0.75, 0.5, 0.25, 0.1;
     Ea   << 1.65e+05, 1.65e+05, 1.65e+05, 1.65e+05, 1.65e+05;
     A    << 1.0e+09, 2.0e+05, 1.0e+09, 1.5e+03, 1.17e+03;
     //A    << 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02, 1.0e+02;
     //beta << 0.0, 1.0, 0.5, 0.25, -0.25;
     //Ea   << 2.5e+05, 2.5e+05, 2.5e+05, 2.5e+05, 2.5e+05;
     //Ea   << 0.0, 0.0, 0.0, 0.0, 0.0;
     //Ea   << 1.75e+05, 1.75e+05, 1.75e+05, 1.75e+05, 1.75e+05;
     
     double kfj; // forward reaction rate coeff.
     double kej; // equilibrium constant
     double kbj; // backward reaction rate coeff.
     double rfj; // forward progress rate
     double rbj; // backware progress rate

     VectorXd rj(Mc); // progress rate

     // Set up s/R - h_RT for each species
     std::vector<double> delta_k(Nc, 0.0);
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
         delta_k[rxn.ProblemInfo.n_species + k - 2] = s_prime[k] / R_universal - h_prime[k] / RT;
     }

     for (int j = 0; j < Mc; j++)
     {
         kfj = A(j) * pow(temperature, beta(j)) * exp(-Ea(j) / RT);
         exp_arg_j = 0.0;
         for (int k = 0; k < Nc; k++)
         {
             exp_arg_j  += nukj(k,j) * delta_k[k];
         }
         kej = pow(pa_RT, gamma(j)) * exp(exp_arg_j);
         kbj = kfj / kej;
         rfj = 1.0;
         rbj = 1.0;
         // Calculate reaction rates
         for (int k = 0; k < Nc; k++)
         {
             rfj *= pow(Yinad[k], nukj_r(k,j));
             rbj *= pow(Yinad[k], nukj_p(k,j));
         }
         rj(j) = kfj * rfj - kbj * rbj;
     }

    return rj;
}


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
  unsigned int n_species = rxn.ProblemInfo.n_species + rxn.ProblemInfo.n_extra;
  unsigned int n_atoms   = rxn.ProblemInfo.n_atoms;
  int dim = n_species + 1;

  if (rxn.ProblemInfo.include_inad)
  {
     dim += n_atoms;
  }

  // Solution vector (units are moles)
  std::vector<double> molar_densities(n_species,0);
  for (unsigned int i = 0; i < n_species; i++)
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
  //std::cout << temperature << std::endl;

  // Get some thermodynamic data for each species
  std::vector<double> h_RT_minus_s_R(n_species);
  rxn.Thermo->h_RT_minus_s_R(temp_cache, h_RT_minus_s_R);

  //std::cout << h_RT_minus_s_R << std::endl;

  // Perform kinetics calculations 
  // (basically just gives RHS of ODEs for each species)
  // (Note that temperature is passed in but the RHS for the
  //  temperature equation is NOT computed by Antioch.  Just
  //  need temperature in reaction rate calculations.)
  std::vector<double> mole_sources(n_species,0.0);
  rxn.Kinetics->compute_mole_sources(
      temperature,     // temperature needed for calculations
      molar_densities, // current concentrations (in moles)
      h_RT_minus_s_R,  // exponent in equilibrium constant for reverse rates
      mole_sources);   // RHS

  /*
  std::cout << molar_densities << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << h_RT_minus_s_R << std::endl;
  std::cout << "\n" << std::endl;
  std::cout << mole_sources << std::endl;
  std::cout << "\n" << std::endl;
  */

  // Set up RHS of each ODE for the species
  for (unsigned int i = 0; i < n_species; i++)
  { // Loop over species
      dYdt[i] = mole_sources[i];
  }

/*
  for (unsigned int k = 0; k < dim; k++)
  {
      std::cout << "dY_" << k << "/dt = " << dYdt[k] << "\n";
  }
  std::cout << "\n";
*/

  // Needed for energy equation
  double Qnum = 0.0;
  double Qden = 0.0;

  double R_universal = Antioch::Constants::R_universal<double>();

  if (rxn.ProblemInfo.include_inad)
  {
     // Create a new vector containing only relevant species
     // This will not be necessary if the user supplies only 
     // the participating species
      std::vector<double> Yinad(rxn.ProblemInfo.n_species, 0.0);
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

     //std::cout << Yinad << std::endl;

     int Mc = 5; // FIXME  Remove hard-coding.

     std::vector<double> h_prime(n_atoms, 0.0);  // Enthalpy for catchalls 
     std::vector<double> cp_prime(n_atoms, 0.0); // Specific heat (constant cp) for catchalls 
     std::vector<double> s_prime(n_atoms, 0.0);  // Entropy for catchalls 

     // Coeffs for catchall thermo
     std::vector<double> alpha_0(n_atoms, 0.0);
     std::vector<double> beta_0 (n_atoms, 0.0);
     std::vector<double> alpha_1(n_atoms, 0.0);
     std::vector<double> alpha_2(n_atoms, 0.0);

     double sig = 0.10;

     // Mean Values
     alpha_0[0] = -3.90372558e+04; //100.5;
     alpha_0[1] = -3.90372558e+04; //150.19;

     // \pm 3std Values
     alpha_0[0] += 3.0 * sig * alpha_0[0];
     alpha_0[1] += 3.0 * sig * alpha_0[1];

     // Mean Values
     beta_0[0] = -17.0857829376; //-0.382;
     beta_0[1] = -17.0857829376; // 24.61;

     // \pm 3std Values
     beta_0[0] += 3.0 * sig * beta_0[0];
     beta_0[1] += 3.0 * sig * beta_0[1];

     // ALPHA_{i1}
     alpha_1[0] = 13.6559654; //1.0e-05;
     alpha_1[1] = 13.6559654; //1.0e-05;

     // \pm 3std Values
     alpha_1[0] += 3.0 * sig * alpha_1[0];
     alpha_1[1] += 3.0 * sig * alpha_1[1];

     // ALPHA_{i2}
     alpha_2[0] = 1.20459536e-03; //1.00e-03;
     alpha_2[1] = 1.20459536e-03; //6.82e-03;

     // \pm 3std Values
     alpha_2[0] += 3.0 * sig * alpha_2[0];
     alpha_2[1] += 3.0 * sig * alpha_2[1];

     for (int m = 0; m < n_atoms; m++)
     {
         h_prime[m]  = alpha_0[m] + 
                       alpha_1[m] * temperature + 
                       alpha_2[m] * temperature * temperature;
         cp_prime[m] = alpha_1[m] + 2.0 * alpha_2[m] * temperature;
         s_prime[m]  = beta_0[m]  + 
                       alpha_1[m] * log(temperature) + 
                       2.0 * alpha_2[m] * temperature;
     }

     Eigen::VectorXd rj(Mc);
     rj = reaction_rate(Yinad, temperature, h_prime, s_prime, rxn);

     // Finally, compute the RHS for these species.
     VectorXd omega_dot_inad;
     MatrixXd nukj   = rxn.ProblemInfo.nukj;
     omega_dot_inad = nukj * rj;
     //std::cout << omega_dot_inad << std::endl;
     //exit(0);

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
         Qnum += h_prime[k] * dYdt[(dim - 1) - n_atoms + k];
         Qden += cp_prime[k] * Y[(dim - 1) - n_atoms + k];
     }
     //std::cout << "Qnum Catchalls = " << Qnum << std::endl;
     //std::cout << "Qnum = " << Qnum << "   " << "Qden = " << Qden << std::endl << std::endl;
     //std::cout << omega_dot_inad << std::endl;

     //std::cout << "omega_dot_inad = " << omega_dot_inad << "\n" << std::endl;
  }

/*
  for (int i = 0; i < dim; i++)
  {
      std::cout << "dx_" << i << "dt = " << dYdt[i] << std::endl;
  }
  std::cout << "\n" << std::endl;
*/

  // Energy equation computations
  std::vector<double> h(n_species, 0.);  // Enthalpy for each species
  std::vector<double> ent(n_species, 0.);  // Entropy for each species
  std::vector<double> cp(n_species, 0.); // Specific pressure for each species

  //double Qnum_real = 0.0;
  for (unsigned int s = 0; s < n_species; s++)
  { // Get numerator and denominator in energy equation (sum of species)
      // Get enthalpy and convert to molar from mass basis
      // Note that R_universal/Rs = Ws where Ws is the molecular weight of species s
      h[s] = rxn.Thermo->h(temp_cache,s) * R_universal / rxn.Chem_mixture->R(s);
      //ent[s] = rxn.Thermo->s_over_R(temp_cache, s) * R_universal;
/*
      std::cout << rxn.Chem_mixture->R(s) << "   " 
                << Antioch::Constants::R_universal<double>() << "   "  
                << Antioch::Constants::R_universal<double>() / rxn.Chem_mixture->R(s) << "   "
                << rxn.Thermo->h(temp_cache,s) << "   " 
                << h[s] <<  std::endl;
*/
      //std::cout << "h_" << s << " = " << h[s] << "   " << "s_" << s << " = " << ent[s] << "   ";
      // get cp and convert from mass to molar
      cp[s] = rxn.Thermo->cp(temp_cache,s) * R_universal / rxn.Chem_mixture->R(s);
      // Numerator in energy equation: h_s * dx_s/dt
      Qnum  += h[s] * mole_sources[s];
      //Qnum_real += h[s] * mole_sources[s];
      // Denominator in energy equation: cp_s * x_s
      Qden += cp[s] * molar_densities[s];
  }

  //std::cout <<  "\n";
  //std::cout << "\n\n\n" << std::endl;
  //std::cout << "omega_dot = " << mole_sources << std::endl;
  //std::cout << "x         = " << molar_densities << std::endl;
  //std::cout << "h         = " << h << std::endl;
  //std::cout << "s         = " << ent << std::endl;
  //std::cout << "cp        = " << cp << std::endl;
  //std::cout << "x = " << molar_densities << std::endl;
  //std::cout << "dxdt = " << mole_sources << std::endl;
  //std::cout << "Qnum total = " << Qnum << std::endl;
  //std::cout << "Qnum real = " << Qnum_real << std::endl << std::endl;
  //std::cout << "Qnum = " << Qnum << "   " << "Qden = " << Qden << std::endl;
  //std::cout << "\n\n\n" << std::endl;
  //exit(0);

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
  double err_abs = 1.0e-12;  // Absolute error tolerance for adaptive stepping
  double err_rel = 1.0e-06;   // Relative error tolerance for adaptive stepping
  gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rkf45, dt, err_abs, err_rel);   

  // Initialize solution field
  double Y[dim];
  for (unsigned int i = 0; i < dim; ++i)
  {
      Y[i] = initialValues[i];
  };

  // Prepare for time integrator to get ignition data
  double t = 0.0; // Initial time
  double dt_samp = 1e-6;
  int Nt = 100000;
  double finalTime;

  bool ignition = false;

  double Tig = 0.0;
  double time_ig = 0.0;

  // Time integration
  for (int i = 0; i < Nt; i++)
  { // Integration loop
      finalTime = i * dt_samp; // Integration time
      while (t < finalTime)
      {
         // std::cout << finalTime << std::endl;
         // Call GSL ODE integrator 
         // t is the current time
         // finalTime is the integration time
         // Y is the solution at time t
         std::cout << "Time = " << t << std::endl;
         gsl_odeiv2_driver_apply( d, &t, finalTime, Y );
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

  std::cout << "\n\n" << "Ignition!" << "\n\n" << std::endl;

  // Reinitialize solution field
  for (int i = 0; i < dim; i++)
  {
      Y[i] = initialValues[i];
  }

  // Prepare for time integrator
  t = 0.0; // Initial time

  // Time integration
  int num_times = timePoints.size();

  int n_atoms = rxn->ProblemInfo.n_atoms;

  std::ofstream progress_rates;
  progress_rates.open("rj.txt");

  std::ofstream reaction_rates;
  reaction_rates.open("omega_k.txt");

  for (unsigned int i = 0; i < num_times; i++)
  { // Integration loop
      finalTime = timePoints[i]; // Integration time
      // Call GSL ODE integrator 
      // t is the current time
      // finalTime is the integration time
      // Y is the solution at time t
      std::cout << "Time = " << t << std::endl;
      int status = gsl_odeiv2_driver_apply( d, &t, finalTime, Y );

      double temperature = Y[dim-1];

      // Create a new vector containing only relevant species
      std::vector<double> Yinad(rxn->ProblemInfo.n_species, 0.0);
      // First just copy participating species to Yinad
      // Skip H and O
      for (int k = 0; k < 2; k++)
      {
          Yinad[k] = Y[k];
      }
      // Skipped H and O so decrement index by 2
      for (int k = 4; k < rxn->ProblemInfo.n_species; k++)
      {
          Yinad[k - 2] = Y[k];
      }
      // Next copy catchall species to Yinad
      for (int k = 0; k < n_atoms; k++)
      {
          Yinad[rxn->ProblemInfo.n_species + k - 2] = Y[(dim - 1) - n_atoms + k];
      }

      int Mc = 5; // FIXME  Remove hard-coding.

      std::vector<double> h_prime(n_atoms, 0.0);  // Enthalpy for catchalls 
      std::vector<double> cp_prime(n_atoms, 0.0); // Specific heat (constant cp) for catchalls 
      std::vector<double> s_prime(n_atoms, 0.0);  // Entropy for catchalls 

      // Coeffs for catchall thermo
      std::vector<double> alpha_0(n_atoms, 0.0);
      std::vector<double> beta_0 (n_atoms, 0.0);
      std::vector<double> alpha_1(n_atoms, 0.0);
      std::vector<double> alpha_2(n_atoms, 0.0);

      double sig = 0.10;

      // Mean Values
      alpha_0[0] = -3.90372558e+04; //100.5;
      alpha_0[1] = -3.90372558e+04; //150.19;

      // \pm 3std Values
      alpha_0[0] += 3.0 * sig * alpha_0[0];
      alpha_0[1] += 3.0 * sig * alpha_0[1];

      // Mean Values
      beta_0[0] = -17.0857829376; //-0.382;
      beta_0[1] = -17.0857829376; // 24.61;

      // \pm 3std Values
      beta_0[0] += 3.0 * sig * beta_0[0];
      beta_0[1] += 3.0 * sig * beta_0[1];

      // ALPHA_{i1}
      alpha_1[0] = 13.6559654; //1.0e-05;
      alpha_1[1] = 13.6559654; //1.0e-05;

      // \pm 3std Values
      alpha_1[0] += 3.0 * sig * alpha_1[0];
      alpha_1[1] += 3.0 * sig * alpha_1[1];

      // ALPHA_{i2}
      alpha_2[0] = 1.20459536e-03; //1.00e-03;
      alpha_2[1] = 1.20459536e-03; //6.82e-03;

      // \pm 3std Values
      alpha_2[0] += 3.0 * sig * alpha_2[0];
      alpha_2[1] += 3.0 * sig * alpha_2[1];

      for (int m = 0; m < n_atoms; m++)
      {
          h_prime[m]  = alpha_0[m] + 
                        alpha_1[m] * temperature + 
                        alpha_2[m] * temperature * temperature;
          cp_prime[m] = alpha_1[m] + 2.0 * alpha_2[m] * temperature;
          s_prime[m]  = beta_0[m]  + 
                        alpha_1[m] * log(temperature) + 
                        2.0 * alpha_2[m] * temperature;
      }

      Eigen::VectorXd rj(Mc);
      VectorXd omega_dot_inad;
      rj = reaction_rate(Yinad, temperature, h_prime, s_prime, *rxn);

      MatrixXd nukj  = rxn->ProblemInfo.nukj;
      omega_dot_inad = nukj * rj;

      progress_rates << t << "   ";
      for (int j = 0; j < Mc; j++)
      {
          progress_rates << rj(j) << "   ";
      }
      progress_rates << "\n";

      reaction_rates << t << "   ";
      for (int k = 0; k < rxn->ProblemInfo.n_species; k++)
      {
          reaction_rates << omega_dot_inad(k) << "   ";
      }
      reaction_rates << "\n";

      // Store results in return field
      for (unsigned int j = 0; j < dim; j++)
      {
          returnValues[dim * i + j] = Y[j];
      }

  } // end loop over time points
  returnValues[dim * num_times] = time_ig;
  returnValues[dim * num_times + 1] = Tig;

  progress_rates.close();
  reaction_rates.close();

  // deallocate memory   
  gsl_odeiv2_driver_free( d );

} // end hydrogenComputeModel
