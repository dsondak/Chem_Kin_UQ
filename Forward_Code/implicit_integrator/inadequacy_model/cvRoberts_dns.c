/*
 * -----------------------------------------------------------------
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>
#include <limits>

// user defined functions
#include "reaction_info.h"
#include "write_data.h"

/* Header files with a description of contents used */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

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

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) DENSE_ELEM(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */
#define RTOL  RCONST(1.0e-08)  /* scalar relative tolerance            */
#define ATOL  RCONST(1.0e-13)  /* scalar relative tolerance            */
#define TIME0    RCONST(0.0)      /* initial time           */
#define DTOUT RCONST(0.000001) /* output time factor     */
#define NOUT  50000            /* number of output times */


/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, const char *funcname, int opt);


/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  realtype reltol, t;
  N_Vector y, abstol;
  void *cvode_mem;
  int flag;

  y = abstol = NULL;
  cvode_mem = NULL;

  const unsigned int n_species        = 7;
  const unsigned int n_inert          = 1;
  const unsigned int n_inad           = 2;
  const unsigned int n_reactions      = 5;
  const unsigned int n_species_tot    = n_species + n_inert + n_inad;
  double             Q                = 5.0e+06; // Heating rate

  const unsigned int n_eq = n_species_tot + 1; // + 1 for temperature

  bool write_to_h5 = true;

  //std::cout << std::numeric_limits<double>::epsilon() << std::endl;
  //exit(0);

  // Now read in the time points
  std::ifstream time_file("time_points.txt");
  std::istream_iterator<double> start(time_file), end;
  std::vector<double> time_points(start, end);
  if (time_file.is_open()) {
     std::istream_iterator<double> start(time_file), end;
     std::vector<double> time_points(start, end);
  }
  else {
     std::cout << "time_points.txt could not be opened." << std::endl;
  }

  if (write_to_h5) {
     create_file("inad_solution.h5", NOUT + 1, n_eq, 1); // for writing out data
  }
  create_file("inad_10pts.h5", time_points.size(), n_eq, 1); // for writing out data


  /* Create serial vector of length n_eq for I.C. and abstol */
  y = N_VNew_Serial(n_eq);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(n_eq); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

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
  species_str_list.push_back("N2");
  species_str_list.push_back("Hp");
  species_str_list.push_back("Op");

  // Get chemistry for species involved in this reaction
  Antioch::ChemicalMixture<double> chem_mixture( species_str_list );

  // Get thermodynamic data
  Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > nasa_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data( nasa_mixture, "nasa7_thermo_inad.xml", Antioch::XML );

  // Prepare for chemical reactions
  Antioch::ReactionSet    <double> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data_xml<double>( "inad_rxn.xml", true, reaction_set );

  // Set up reactions and thermodynamics 
  Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0); // Reactions
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > thermo(nasa_mixture); // Thermodynamics

  reaction_info rxnMain(&chem_mixture, &reaction_set, &thermo, &nasa_mixture, &kinetics, n_species, n_inert, n_inad, 
                        n_eq, n_reactions, Q);

  /* Initialize y */
  std::vector<double> concs(n_species_tot, std::numeric_limits<double>::epsilon());
  concs[0] = 2.0;
  concs[1] = 1.0;
  concs[7] = 3.78;
  for (unsigned int i = 0; i < n_species_tot; i++) {
      Ith(y,i+1) = concs[i];
  }
  Ith(y,n_eq) = 450.0; // temperature

  /* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
  for (unsigned int i = 0; i < n_eq; i++) {
      Ith(abstol,i+1) = ATOL;
  }

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  
  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time TIME0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, TIME0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

  /* Set stability limit detection algorithm */
  //flag = CVodeSetStabLimDet(cvode_mem, true); // Default is false
  //if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

  /* Set minimum time step size */
  //realtype hmin = 1.0e-16;
  //flag = CVodeSetMinStep(cvode_mem, hmin);
  //if (check_flag(&flag, "CVodeSetMinStep", 1)) return(1);

  /* Set maximum time step size */
  //realtype hmax = 5.0e-07; // 0.0 means \infty
  //realtype hmax = 1.0e-07;
  //flag = CVodeSetMaxStep(cvode_mem, hmax);
  //if (check_flag(&flag, "CVodeSetMaxStep", 1)) return(1);

  /* Set maximum number of nonlinear solver iterations */
  //flag = CVodeSetMaxNonlinIters(cvode_mem, 100);
  //if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

  /* Reset nonlinear convergence test safety factor */
  //flag = CVodeSetNonlinConvCoef(cvode_mem, 0.1);
  //if (check_flag(&flag, "CvodeSetNonlinConvCoef", 1)) return(1);

  /* Set max number of nonlinear solver convergence failures */
  //flag = CVodeSetMaxConvFails(cvode_mem, 1000); // default = 10
  //if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

  /* Set max number of error test failures */
  //flag = CVodeSetMaxErrTestFails(cvode_mem, 1000000); // default = 7
  //if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return(1);

  /* Set maximum order of the linear multistep method */
  flag = CVodeSetMaxOrd(cvode_mem, 1);
  if (check_flag(&flag, "CVodeSetMaxOrd", 1)) return(1);

  /* Set maximum number of internal steps */
  flag = CVodeSetMaxNumSteps(cvode_mem, 1000000000);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, n_eq);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

  /* Set user data */
  flag = CVodeSetUserData(cvode_mem, &rxnMain);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  //printf(" \n5-species kinetics problem\n\n");
  std::cout << n_species_tot << " kinetics problem." << std::endl;

  // Implement constraints to reset the entropy coefficients
  double R_universal = Antioch::Constants::R_universal<double>();
  double molecular_weight;

  std::vector<double> inad_constraints(n_inad, 0.0);
  std::vector<double> alphas(3*n_inad, 0.0);

  inad_constraints[0] = 160.992;
  inad_constraints[1] = 162.901;

  alphas[0] = -5.93184;
  alphas[1] =  2.45189;
  alphas[2] =  2.08591e-04;
  alphas[3] = -4.15409;
  alphas[4] =  2.44766;
  alphas[5] =  4.89017e-04;

  std::vector<double> h_prime_Tmin(n_inad, 0.0);
  std::vector<double> h_prime_Tmax(n_inad, 0.0);
  std::vector<double> sTmin(n_inad, 0.0);
  std::vector<double> sTmax(n_inad, 0.0);

  double Tmin   = 450.0; // Temperature better not go below this
  double Tmax   = 3500.0; // Temperature will hopefully not exceed this
  double RT_min = R_universal * Tmin;
  double RT_max = R_universal * Tmax;

  // Set min and max temperature caches
  Antioch::TempCache<double> temp_cache_Tmin(Tmin);
  Antioch::TempCache<double> temp_cache_Tmax(Tmax);

  // Get some thermodynamic data for each species
  std::vector<double> h_RT_minus_s_R_Tmin(n_species_tot, 0.0);
  std::vector<double> h_RT_minus_s_R_Tmax(n_species_tot, 0.0);
  rxnMain.Thermo->h_RT_minus_s_R(temp_cache_Tmin, h_RT_minus_s_R_Tmin);
  rxnMain.Thermo->h_RT_minus_s_R(temp_cache_Tmax, h_RT_minus_s_R_Tmax);

  // Step 1:  Compute enthalpy and specific heat for virtual species
  for (unsigned int k = 0; k < n_inad; k++) {
     molecular_weight = R_universal / rxnMain.Chem_mixture->R(n_species + n_inert + k);
     h_prime_Tmin[k]  = rxnMain.Thermo->h(temp_cache_Tmin,  n_species + n_inert + k) * molecular_weight;
     h_prime_Tmax[k]  = rxnMain.Thermo->h(temp_cache_Tmax,  n_species + n_inert + k) * molecular_weight;
  }
  // Step 3:  Compute temperature-dependent terms in entropy
  for (unsigned int k = 0; k < n_inad; k++) {
     sTmin[k] = alphas[3*k + 1] * log(Tmin) + 2.0 * alphas[3*k + 2] * Tmin;
     sTmax[k] = alphas[3*k + 1] * log(Tmax) + 2.0 * alphas[3*k + 2] * Tmax;
  }
  // Step 4:  Calculate first constant term in entropy
  double c1_Tmin = -sTmin[0] + (h_prime_Tmin[0] - 0.5 * h_RT_minus_s_R_Tmin[0] * RT_min) / Tmin;
  double c1_Tmax = -sTmax[0] + (h_prime_Tmax[0] - 0.5 * h_RT_minus_s_R_Tmax[0] * RT_max) / Tmax;

  // Find the minimum of c1
  double c1;
  if (c1_Tmax < c1_Tmin) {
     c1 = c1_Tmax;
  }
  else if (c1_Tmin < c1_Tmax) {
     c1 = c1_Tmin;
  }
  else {
     c1 = c1_Tmin;
  }

  std::vector<double> betas(n_inad, 0.0);

  // Get the first entropy coefficient
  betas[0] = (c1 - inad_constraints[0]) / R_universal;

  // Now we have to go back and do it again since the remaining 
  // constraints depend on beta_0
  std::vector<double> constraints(4, 0.0); // Vector of constraints

  // Compute c3 and c4 which also have their minima at Tmin
  constraints[1] = -betas[0] - sTmin[0] - sTmin[1] + 
          (h_prime_Tmin[0] + h_prime_Tmin[1] - h_RT_minus_s_R_Tmin[4] * RT_min) / Tmin;
  constraints[2] = 0.5 * (-betas[0] - sTmin[0] - 2.0 * sTmin[1] + 
                (h_prime_Tmin[0] + 2.0 * h_prime_Tmin[1] - h_RT_minus_s_R_Tmin[5] * RT_min) / Tmin);

  constraints[0] = -sTmax[1] + (h_prime_Tmax[1] - 0.5 * h_RT_minus_s_R_Tmax[1] * RT_max) / Tmax;
  constraints[3] = -2.0 * betas[0] - 2.0 * sTmax[0] - sTmax[1] + 
          (2.0 * h_prime_Tmax[0] + h_prime_Tmax[1] - h_RT_minus_s_R_Tmax[6] * RT_max) / Tmax;

  // Step 6:  Find the most restrictive constraint
  double min_constraint = *std::min_element(constraints.begin(), constraints.end());
  // Step 7:  Calculate the second constant term in entropy
  betas[1] = (min_constraint - inad_constraints[1]) / R_universal;
  // Step 8:  Reset the appropriate coefficients in the NASA polynomial
  for (unsigned int k = 0; k < n_inad; k++) {
     rxnMain.NASAMixture->set_curve_fit_coefficient(n_species + n_inert + k, 0, 6, betas[k]);
     rxnMain.NASAMixture->set_curve_fit_coefficient(n_species + n_inert + k, 1, 6, betas[k]);
  }

  /* Set up time and solution vectors */
  std::vector<double> timef(NOUT + 1, 0.0);
  std::vector<double> solutionf(n_eq * (NOUT + 1) + 2, 0.0);

  /* Initialize time and solution vectors */
  timef[0] = TIME0;
  for (unsigned int j = 0; j < n_eq; j++) {
      solutionf[j] = Ith(y,j+1);
  }

  /* Get ignition data */
  double time_ig = 0.0;
  double T_ig = 0.0;
  bool no_ig = true;

  for (unsigned int iout = 1; iout <= NOUT; iout++) {

    realtype tout = iout * DTOUT;

    try {
       flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    }
    catch(...){
       if (write_to_h5) {
          write_file(timef, solutionf, "inad_solution.h5", NOUT + 1, n_eq, 0, 1.0, 5.0e+06);
       }
       else {
          std::ofstream txt_solution;
          txt_solution.open("inad_solution.txt");
          for (unsigned int i = 0; i <= NOUT; i++) {
              txt_solution << timef[i] << "   ";
              for (unsigned int j = 0; j < n_eq; j++) {
                  txt_solution << solutionf[i * n_eq + j] << "   ";
              }
              txt_solution << "\n";
          }
          txt_solution.close();
       }
       exit(0);
    }

    std::cout << "t = " << t << std::endl;

    if (check_flag(&flag, "CVode", 1)) break;

    if (no_ig) {
       if (Ith(y,n_eq) >= 1500.0) {
          time_ig = t;
          T_ig = Ith(y,n_eq);
          no_ig = false;
       }
    }

    /* Store time and solution */
    timef[iout] = tout;

    for (unsigned int j = 0; j < n_eq; j++) {
        solutionf[iout * n_eq + j] = Ith(y,j+1);
    }

  }

  /* Ignition data */
  solutionf[n_eq * (NOUT + 1)]     = time_ig;
  solutionf[n_eq * (NOUT + 1) + 1] = T_ig;

  /* Write solution to file */
  if (write_to_h5) {
     write_file(timef, solutionf, "inad_solution.h5", NOUT + 1, n_eq, 0, 1.0, 5.0e+06);
  }
  else {
     std::ofstream txt_solution;
     txt_solution.open("inad_solution.txt");
     for (unsigned int i = 0; i <= NOUT; i++) {
         txt_solution << timef[i] << "   ";
         for (unsigned int j = 0; j < n_eq; j++) {
             txt_solution << solutionf[i * n_eq + j] << "   ";
         }
         txt_solution << "\n";
     }
     txt_solution.close();
  }

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// Now just evaluate the solution at a few pre-determined times
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  /* Re-Initialize y */
  for (unsigned int i = 0; i < n_species_tot; i++) {
      Ith(y,i+1) = concs[i];
  }
  Ith(y,n_eq) = 450.0; // temperature

  flag = CVodeReInit(cvode_mem, TIME0, y);
  if (check_flag(&flag, "CVodeReInit", 1)) return(1);

  /* Set up time and solution vectors */
  std::vector<double> time(time_points.size(), 0.0);
  std::vector<double> solution(n_eq * time_points.size() + 2, 0.0);

  double time_ig_truth = 0.019616;
  double dt_ig = time_ig_truth - time_ig;
  std::cout << "dt_ig = " << dt_ig << std::endl;
  for (unsigned int iout = 0; iout < time_points.size(); iout++) {

    realtype tout = time_points[iout] - dt_ig;

    std::cout << "t = " << tout << std::endl;

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    if (check_flag(&flag, "CVode", 1)) break;

    /* Store time and solution */
    time[iout] = tout;

    for (unsigned int j = 0; j < n_eq; j++) {
        solution[iout * n_eq + j] = Ith(y,j+1);
    }

  }

  /* Ignition data */
  solution[n_eq * time_points.size()] = time_ig;
  solution[n_eq * time_points.size() + 1] = T_ig;

  /* Write solution to file */
  write_file(time, solution, "inad_10pts.h5", time_points.size(), n_eq, 0, 2.0, 5.0e+06);

  /* Print some final statistics */
  PrintFinalStats(cvode_mem);

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */

/*
 * f routine. Compute function f(t,y). 
 */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{

  reaction_info rxn = *(reaction_info *) user_data;

  const unsigned int n_species     = rxn.ProblemInfo.n_species;
  const unsigned int n_inert       = rxn.ProblemInfo.n_inert;
  const unsigned int n_inad        = rxn.ProblemInfo.n_inad;
  const unsigned int n_eq          = rxn.ProblemInfo.n_eq;
  const unsigned int n_species_tot = n_species + n_inert + n_inad;

  // Solution vector (units are moles)
  std::vector<double> molar_densities(n_species_tot,0);
  for (unsigned int i = 0; i < n_species_tot; i++)
  { // Populate molar_densities vector
      molar_densities[i] = Ith(y,i+1);
      if (molar_densities[i] < 0.0) {
         molar_densities[i] = 0.0;
      }
  }

  // Prepare Antioch to get chemistry and thermodynamic information
  double temperature = Ith(y,n_eq);
  Antioch::TempCache<double> temp_cache(temperature);

  // Get some thermodynamic data for each species
  std::vector<double> h_RT_minus_s_R(n_species_tot, 0.0);
  rxn.Thermo->h_RT_minus_s_R(temp_cache, h_RT_minus_s_R);

  // Perform kinetics calculations 
  std::vector<double> mole_sources(n_species_tot,0.0);
  rxn.Kinetics->compute_mole_sources(
      temperature,     // temperature needed for calculations
      molar_densities, // current concentrations (in moles)
      h_RT_minus_s_R,  // exponent in equilibrium constant for reverse rates
      mole_sources);   // RHS

  for (unsigned int i = 0; i < n_species_tot; i++) { 
      // Populate molar_densities vector
      Ith(ydot,i+1) = mole_sources[i];
  }

  // Energy equation
  double h_dot_mole_sources = 0.0;
  double cp_dot_species     = 0.0;
  double molecular_weight   = 0.0;
  double R_universal        = Antioch::Constants::R_universal<double>();

  for (unsigned int k = 0; k < n_species_tot; k++) {
      molecular_weight    = R_universal / rxn.Chem_mixture->R(k);
      h_dot_mole_sources += rxn.Thermo->h(temp_cache, k) * molecular_weight * mole_sources[k];
      cp_dot_species     += rxn.Thermo->cp(temp_cache, k) * molecular_weight * molar_densities[k];
  }

  Ith(ydot,n_eq) = (-h_dot_mole_sources + rxn.ProblemInfo.Q) / cp_dot_species;

  int return_val = 0;

  return(return_val);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

  reaction_info rxn = *(reaction_info *) user_data;

  const unsigned int n_species     = rxn.ProblemInfo.n_species;
  const unsigned int n_inert       = rxn.ProblemInfo.n_inert;
  const unsigned int n_inad        = rxn.ProblemInfo.n_inad;
  const unsigned int n_eq          = rxn.ProblemInfo.n_eq;
  const unsigned int n_species_tot = n_species + n_inert + n_inad;

  // Solution vector (units are moles)
  std::vector<double> molar_densities(n_species_tot,0);
  for (unsigned int i = 0; i < n_species_tot; i++)
  { // Populate molar_densities vector
      molar_densities[i] = Ith(y,i+1);
      if (molar_densities[i] < 0.0) {
         molar_densities[i] = 0.0;
      }
  }

  // Prepare Antioch to get chemistry and thermodynamic information
  double temperature = Ith(y,n_eq);
  Antioch::TempCache<double> temp_cache(temperature);

  // Gibbs free energy for each species
  std::vector<double> h_RT_minus_s_R(n_species_tot, 0.0);
  rxn.Thermo->h_RT_minus_s_R(temp_cache, h_RT_minus_s_R);

  // Derivative of Gibbs free energy for each species wrt T
  std::vector<double> dh_RT_minus_s_R_dT(n_species_tot, 0.0);
  rxn.Thermo->dh_RT_minus_s_R_dT(temp_cache, dh_RT_minus_s_R_dT);

  // Derivatives of RHS wrt T and species
  std::vector<double> dmole_sources_dT(n_species_tot, 0.0);
  std::vector<std::vector<double> > dmole_sources_dX_s(n_species_tot);

  for (unsigned int k = 0; k < n_species_tot; k++) {
      dmole_sources_dX_s[k].resize(n_species_tot);
  }

  // Perform kinetics calculations 
  std::vector<double> mole_sources(n_species_tot, 0.0);
  rxn.Kinetics->compute_mole_sources_and_derivs(
      temperature,     // temperature needed for calculations
      molar_densities, // current concentrations (in moles)
      h_RT_minus_s_R,  // exponent in equilibrium constant for reverse rates
      dh_RT_minus_s_R_dT, // derivative of Keq. exp.
      mole_sources,    // RHS
      dmole_sources_dT, // Derivative of RHS wrt temperature
      dmole_sources_dX_s // Derivative of RHS wrt species
  );

  for (unsigned int i = 0; i < n_species_tot; i++) {
      for (unsigned int j = 0; j < n_species_tot; j++) {
          IJth(J,i+1,j+1) = dmole_sources_dX_s[i][j];
      }
      IJth(J,i+1,n_eq) = dmole_sources_dT[i];
  }

  // Now do Jacobian for energy equation
  double h_dot_mole_sources = 0.0;
  double cp_dot_mole_sources= 0.0;
  double cp_dot_species     = 0.0;
  double h_dot_J_T          = 0.0;
  double dcp_dT_dot_species = 0.0;
  double molecular_weight   = 0.0;
  double R_universal        = Antioch::Constants::R_universal<double>();

  for (unsigned int k = 0; k < n_species_tot; k++) {
      molecular_weight    = R_universal / rxn.Chem_mixture->R(k);
      h_dot_mole_sources += rxn.Thermo->h(temp_cache, k) * molecular_weight * mole_sources[k];
      cp_dot_mole_sources+= rxn.Thermo->cp(temp_cache, k) * molecular_weight * mole_sources[k];
      cp_dot_species     += rxn.Thermo->cp(temp_cache, k) * molecular_weight * molar_densities[k];
      h_dot_J_T          += rxn.Thermo->h(temp_cache, k) * molecular_weight * IJth(J,k+1,n_eq);
      dcp_dT_dot_species += rxn.Thermo->dcp_dT(temp_cache, k) * molecular_weight * molar_densities[k];
  }

  double cp_dot_species_2 = cp_dot_species * cp_dot_species;

  // Energy RHS derivatives wrt species
  double num1;
  double num2;
  for (unsigned int j = 0; j < n_species_tot; j++) {
      double h_dot_J = 0.0;
      for (unsigned int k = 0; k < n_species_tot; k++) {
          h_dot_J += rxn.Thermo->h(temp_cache, k) * molecular_weight * IJth(J,k+1,j+1);
      }
      num1 = -cp_dot_species * h_dot_J;
      num2 = (-h_dot_mole_sources + rxn.ProblemInfo.Q) * rxn.Thermo->cp(temp_cache, j);
      IJth(J,n_eq, j+1) = (num1 - num2) / cp_dot_species_2;
  }

  // Energy RHS derivative wrt temperature
  num1 = cp_dot_species * (h_dot_J_T + cp_dot_mole_sources);
  num2 = (-h_dot_mole_sources + rxn.ProblemInfo.Q) * dcp_dT_dot_species;
  IJth(J,n_eq,n_eq) = (num1 - num2) / cp_dot_species_2;

  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}
