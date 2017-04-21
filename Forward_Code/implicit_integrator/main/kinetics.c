/***********************************************************
 * Forward code for chemical kinetics.
 *
 * Solves a system of ODEs
 * 
 *      dx/dt = f(x, T)
 *      dT/dt = fT(x,T)
 * 
 * where x \in R^Ns is a vector of 
 * species with units of moldes and 
 * T \in R is a scalar denoting 
 * temperature in Kelvin.
 *
 * The system is solved with CVODE 
 * using adaptive order and timestep 
 * BDFs methods.
 *
 * The RHS of the species equations 
 * is provided by Antioch.  The RHS
 * for the energy equation is computed 
 * manually.
 *
 * Solutions are written out to HDF5 
 * files.
 *
 * LINES TO CHANGE BEFORE RUNNING:
 * Input parameters block starting at line 144
 *
 *
 * TO DO:  MAKE AN INPUT FILE
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>

// user defined functions
#include "reaction_info.h"
#include "write_data.h"

/* Header files with a description of contents used */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <nlopt.h>                   /* nlopt optimization package */

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
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.000001) /* first output time      */
#define DTOUT RCONST(0.000001) /* output time factor     */
#define NOUT  50000            /* number of output times */


/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int Jac(long int N, realtype t,
               N_Vector y, N_Vector fy, DlsMat J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int mass_balance(unsigned int n_species, unsigned int n_atoms, 
                 std::vector<double> Amat, std::vector<double> bvec, 
                 double *x1);

typedef struct {
    double* x1;
} my_f_data;

double myfunc(unsigned int n, const double *x, double *grad, void *my_func_data);

typedef struct {
    double* b;
    double* A;
} my_constraint_data;

void myconstraint(unsigned int m, double *result, unsigned int n, const double* x, 
                  double* grad, void* f_data);

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

  /*****************************
   *
   *     INPUT PARAMETERS
   *
   ******************************/
  // Input and output filenames
  std::string thermo_fname("nasa7_thermo_reduced.xml"); // thermo data file
  std::string reaction_set_fname("five_rxn.xml");    // reaction set file

  // Ugly C syntax b/c using C HDF5 interface
  //char data_fname[21];
  char data_fname[20];
  strcpy(data_fname, "reduced_solution.h5");  // file to write solution to

  const unsigned int n_species        = 7;  // number of species
  const unsigned int n_inert          = 1;  // e.g. N2
  const unsigned int n_inad           = 0;  // virtual species
  const unsigned int n_atoms          = 2;  // types of atoms in the system
  const unsigned int n_reactions      = 5; // number of reactions
  const unsigned int n_species_tot    = n_species + n_inert + n_inad;
  double             Q                = 5.0e+06; // heating rate

  double init_H2 = 2.0;   // Initial moles of H2
  double init_O2 = 1.0;   // Initial moles of O2
  double init_N2 = 3.78;  // Initial moles of N2
  double init_T  = 450.0; // Initial temperature (in K)

  const unsigned int n_eq = n_species_tot + 1; // + 1 for temperature

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
  //species_str_list.push_back("H2O2"); // remove this for the reduced model
  species_str_list.push_back("N2");

  // Number of atoms of each type in the system
  std::vector<double> bvec(n_atoms + 1, 0.0);
  bvec[0] = 2.0 * init_H2;
  bvec[1] = 2.0 * init_O2;
  bvec[2] = 2.0 * init_N2;

  // Species to atoms transition matrix
  std::vector<double> Amat((n_atoms + 1) * n_species_tot, 0.0);
  // For detailed model
  //Amat[0] = 2.0;
  //Amat[1] = 0.0;
  //Amat[2] = 1.0;
  //Amat[3] = 0.0;
  //Amat[4] = 1.0;
  //Amat[5] = 1.0;
  //Amat[6] = 2.0;
  //Amat[7] = 2.0;
  //Amat[8] = 0.0;

  //Amat[9]  = 0.0;
  //Amat[10] = 2.0;
  //Amat[11] = 0.0;
  //Amat[12] = 1.0;
  //Amat[13] = 1.0;
  //Amat[14] = 2.0;
  //Amat[15] = 1.0;
  //Amat[16] = 2.0;
  //Amat[17] = 0.0;

  //Amat[18] = 0.0;
  //Amat[19] = 0.0;
  //Amat[20] = 0.0;
  //Amat[21] = 0.0;
  //Amat[22] = 0.0;
  //Amat[23] = 0.0;
  //Amat[24] = 0.0;
  //Amat[25] = 0.0;
  //Amat[26] = 2.0;

  // For reduced model
  Amat[0] = 2.0;
  Amat[1] = 0.0;
  Amat[2] = 1.0;
  Amat[3] = 0.0;
  Amat[4] = 1.0;
  Amat[5] = 1.0;
  Amat[6] = 2.0;
  Amat[7] = 0.0;

  Amat[8] = 0.0;
  Amat[9] = 2.0;
  Amat[10] = 0.0;
  Amat[11] = 1.0;
  Amat[12] = 1.0;
  Amat[13] = 2.0;
  Amat[14] = 1.0;
  Amat[15] = 0.0;

  Amat[16] = 0.0;
  Amat[17] = 0.0;
  Amat[18] = 0.0;
  Amat[19] = 0.0;
  Amat[20] = 0.0;
  Amat[21] = 0.0;
  Amat[22] = 0.0;
  Amat[23] = 2.0;

  // Create output file
  bool write_to_h5 = true;

  if (write_to_h5) {
     create_file(data_fname, NOUT + 1, n_eq, 1); // for writing out data
  }

  /* END INPUT PARAMETERS */

  /* Create serial vector of length n_eq for I.C. and abstol */
  y = N_VNew_Serial(n_eq);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(n_eq); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

  // Get chemistry for species involved in this reaction
  Antioch::ChemicalMixture<double> chem_mixture( species_str_list );

  // Get thermodynamic data
  Antioch::NASAThermoMixture<double, Antioch::NASA7CurveFit<double> > nasa_mixture( chem_mixture );
  Antioch::read_nasa_mixture_data( nasa_mixture, thermo_fname, Antioch::XML );

  // Prepare for chemical reactions
  Antioch::ReactionSet    <double> reaction_set( chem_mixture );
  Antioch::read_reaction_set_data_xml<double>( reaction_set_fname, true, reaction_set );

  // Set up reactions and thermodynamics 
  Antioch::KineticsEvaluator<double> kinetics(reaction_set, 0); // Reactions
  Antioch::NASAEvaluator<double, Antioch::NASA7CurveFit<double> > thermo(nasa_mixture); // Thermodynamics

  reaction_info rxnMain(&chem_mixture, &reaction_set, &thermo, &nasa_mixture, &kinetics, Amat, bvec, 
                        n_species, n_inert, n_inad, n_atoms, n_eq, n_reactions, Q);

  /* Initialize y */
  std::vector<double> init_moles(n_species_tot, 0.0);
  init_moles[0] = init_H2;
  init_moles[1] = init_O2;
  init_moles[n_species + n_inad] = init_N2; // N2
  for (unsigned int i = 0; i < n_species_tot; i++) {
      Ith(y,i+1) = init_moles[i];
  }
  Ith(y,n_eq) = init_T; // temperature

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
   * user's right hand side function in y'=f(t,y), the inital time T0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, T0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

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
  printf("%5.1i %s\n", n_species_tot, "species kinetics problem");

  /* Set up time and solution vectors */
  std::vector<double> timef(NOUT + 1, 0.0);
  std::vector<double> solutionf(n_eq * (NOUT + 1) + 2, 0.0);

  /* Initialize time and solution vectors */
  timef[0] = T0;
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
          write_file(timef, solutionf, data_fname, NOUT + 1, n_eq, 0, 1.0, 5.0e+06);
       }
       else {
          std::ofstream txt_solution;
          txt_solution.open("detailed_solution_cvode.txt");
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

    printf("TIME = %25.16e\n", t);

    if (check_flag(&flag, "CVode", 1)) break;

    if (no_ig) {
       if (Ith(y,n_eq) >= 1500.0) {
          time_ig = t;
          T_ig = Ith(y,n_eq);
          no_ig = false;
          //rxnMain.ProblemInfo.Q = 0.0;
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
     write_file(timef, solutionf, data_fname, NOUT + 1, n_eq, 0, 1.0, 5.0e+06);
  }
  else {
     std::ofstream txt_solution;
     txt_solution.open("detailed_solution_cvode.txt");
     for (unsigned int i = 0; i <= NOUT; i++) {
         txt_solution << timef[i] << "   ";
         for (unsigned int j = 0; j < n_eq; j++) {
             txt_solution << solutionf[i * n_eq + j] << "   ";
         }
         txt_solution << "\n";
     }
     txt_solution.close();
  }

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
  const unsigned int n_atoms       = rxn.ProblemInfo.n_atoms;
  const unsigned int n_eq          = rxn.ProblemInfo.n_eq;
  const unsigned int n_species_tot = n_species + n_inert + n_inad;

  // Temporary solution vector (units are moles)
  int neg_flag = 0;
  double temp_moles[n_species_tot];
  for (unsigned int i = 0; i < n_species_tot; i++)
  {
      temp_moles[i] = Ith(y,i+1);
      if (temp_moles[i] < 0.0) {
         neg_flag += 1;
      }
  }

  if (neg_flag != 0) {

     int mass_flag = mass_balance(n_species_tot, n_atoms, rxn.Amat, rxn.bvec, temp_moles);

     if (mass_flag < 0) {
        printf("Mass balance incomplete.  Error Code: %2.1i.\n", mass_flag);
        if (mass_flag == -4) {
           printf("Failure was due to roundoff errors.\n");
           int neg_flag = 0;
           for (unsigned int k = 0; k < n_species_tot; k++) {
               if (temp_moles[k] < 0) {
                  neg_flag += 1;
               }
           }
           if (neg_flag != 0) {
              printf("Mass balance failed.  Here are the species:\n");
              for (unsigned int i = 0; i < n_species_tot; i++)
              {
                  printf("x_%2.2i = %25.16e\n", i+1, temp_moles[i]);
              }
              exit(0);
           }
           else {
               printf("Mass balance is okay.\n");
           }
        }
        else {
           printf("Mass balance failed.  Error Code: %2.1i.\n", mass_flag);
           printf("Here is the species set:");
           for (unsigned int i = 0; i < n_species_tot; i++)
           {
               printf("x_%2.2i = %25.16e\n", i+1, temp_moles[i]);
           }
           exit(0);
        }
     } // end mass_flag check
  } // end neg_flag check

  // Compute volume from ideal gas law
  double total_moles = 0.0;
  for (unsigned int k = 0; k < n_species_tot; k++) {
      total_moles += temp_moles[k];
  }
  double R_universal = Antioch::Constants::R_universal<double>();
  double press = 1.0e+05;
  double temperature = Ith(y,n_eq);
  Antioch::TempCache<double> temp_cache(temperature);
  double V = total_moles * R_universal * temperature / press;
  
  // Solution vector (units are moles)
  std::vector<double> molar_densities(n_species_tot,0);
  for (unsigned int i = 0; i < n_species_tot; i++)
  { // Populate molar_densities vector
      molar_densities[i] = temp_moles[i] / V;
  }

  // Get some thermodynamic data for each species
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
      temperature,        // temperature needed for calculations
      molar_densities,    // current concentrations (in moles)
      h_RT_minus_s_R,     // exponent in equilibrium constant for reverse rates
      dh_RT_minus_s_R_dT, // derivative of Keq. exp.
      mole_sources,       // RHS
      dmole_sources_dT,   // Derivative of RHS wrt temperature
      dmole_sources_dX_s  // Derivative of RHS wrt species
  );

  // RHS for species molar concentrations
  for (unsigned int k = 0; k < n_species_tot; k++) {
      Ith(ydot, k+1) = mole_sources[k] * V;
  }

  // Energy equation
  double h_dot_mole_sources = 0.0;
  double cp_dot_species     = 0.0;
  double molecular_weight   = 0.0;

  for (unsigned int k = 0; k < n_species_tot; k++) {
      molecular_weight    = R_universal / rxn.Chem_mixture->R(k);
      h_dot_mole_sources += rxn.Thermo->h(temp_cache, k) * molecular_weight * mole_sources[k];
      cp_dot_species     += rxn.Thermo->cp(temp_cache, k) * molecular_weight * Ith(y,k+1);
  }

  Ith(ydot,n_eq) = (-h_dot_mole_sources * V + rxn.ProblemInfo.Q) / cp_dot_species;

  return(0);
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
  const unsigned int n_atoms       = rxn.ProblemInfo.n_atoms;
  const unsigned int n_eq          = rxn.ProblemInfo.n_eq;
  const unsigned int n_species_tot = n_species + n_inert + n_inad;

  // Temporary solution vector (units are moles)
  int neg_flag = 0;
  double temp_moles[n_species_tot];
  for (unsigned int i = 0; i < n_species_tot; i++)
  {
      temp_moles[i] = Ith(y,i+1);
      if (temp_moles[i] < 0.0) {
         neg_flag += 1;
      }
  }

  if (neg_flag != 0) {

     int mass_flag = mass_balance(n_species_tot, n_atoms, rxn.Amat, rxn.bvec, temp_moles);

     if (mass_flag < 0) {
        printf("Mass balance incomplete.  Error Code: %2.1i.\n", mass_flag);
        if (mass_flag == -4) {
           printf("Failure was due to roundoff errors.\n");
           int neg_flag = 0;
           for (unsigned int k = 0; k < n_species_tot; k++) {
               if (temp_moles[k] < 0) {
                  neg_flag += 1;
               }
           }
           if (neg_flag != 0) {
              printf("Mass balance failed.  Here are the species:\n");
              for (unsigned int i = 0; i < n_species_tot; i++)
              {
                  printf("x_%2.2i = %25.16e\n", i+1, temp_moles[i]);
              }
              exit(0);
           }
           else {
               printf("Mass balance is okay.\n");
           }
        }
        else {
           printf("Mass balance failed.  Error Code: %2.1i.\n", mass_flag);
           printf("Here is the species set:");
           for (unsigned int i = 0; i < n_species_tot; i++)
           {
               printf("x_%2.2i = %25.16e\n", i+1, temp_moles[i]);
           }
           exit(0);
        }
     }
  }

  // Compute volume from ideal gas law
  double total_moles = 0.0;
  for (unsigned int k = 0; k < n_species_tot; k++) {
      total_moles += temp_moles[k];
  }
  double R_universal = Antioch::Constants::R_universal<double>();
  double press = 1.0e+05;
  double temperature = Ith(y,n_eq);
  Antioch::TempCache<double> temp_cache(temperature);
  double V = total_moles * R_universal * temperature / press;
  
  // Solution vector (units are moles)
  std::vector<double> molar_densities(n_species_tot,0);
  for (unsigned int i = 0; i < n_species_tot; i++)
  { // Populate molar_densities vector
      molar_densities[i] = temp_moles[i] / V;
  }

  // Get derivatives of volume
  double dV_dxj = R_universal * temperature / press;
  double dV_dT  = R_universal * total_moles / press;

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
      temperature,        // temperature needed for calculations
      molar_densities,    // current concentrations (in moles)
      h_RT_minus_s_R,     // exponent in equilibrium constant for reverse rates
      dh_RT_minus_s_R_dT, // derivative of Keq. exp.
      mole_sources,       // RHS
      dmole_sources_dT,   // Derivative of RHS wrt temperature
      dmole_sources_dX_s  // Derivative of RHS wrt species
  );

  for (unsigned int k = 0; k < n_species_tot; k++) {
      for (unsigned int j = 0; j < n_species_tot; j++) {
          IJth(J, k+1, j+1) = dmole_sources_dX_s[k][j] * V + mole_sources[k] * dV_dxj;
      }
      IJth(J, k+1, n_eq) = dmole_sources_dT[k] * V + mole_sources[k] * dV_dT;
  }

  // Now do Jacobian for energy equation
  double h_dot_mole_sources = 0.0;
  double cp_dot_mole_sources= 0.0;
  double cp_dot_species     = 0.0;
  double h_dot_J_T          = 0.0;
  double dcp_dT_dot_species = 0.0;
  double molecular_weight   = 0.0;

  for (unsigned int k = 0; k < n_species_tot; k++) {
      molecular_weight    = R_universal / rxn.Chem_mixture->R(k);
      h_dot_mole_sources += rxn.Thermo->h(temp_cache, k) * molecular_weight * mole_sources[k];
      cp_dot_mole_sources+= rxn.Thermo->cp(temp_cache, k) * molecular_weight * mole_sources[k];
      cp_dot_species     += rxn.Thermo->cp(temp_cache, k) * molecular_weight * molar_densities[k];
      h_dot_J_T          += rxn.Thermo->h(temp_cache, k) * molecular_weight * 
                              (mole_sources[k] * dV_dT + dmole_sources_dT[k] * V) ;
      dcp_dT_dot_species += rxn.Thermo->dcp_dT(temp_cache, k) * molecular_weight * Ith(y, k+1);
  }

  double cp_dot_species_2 = cp_dot_species * cp_dot_species;

  // Energy RHS derivatives wrt species
  double num1;
  double num2;
  for (unsigned int j = 0; j < n_species_tot; j++) {
      double h_dot_J = 0.0;
      for (unsigned int k = 0; k < n_species_tot; k++) {
          molecular_weight    = R_universal / rxn.Chem_mixture->R(k);
          h_dot_J += rxn.Thermo->h(temp_cache, k) * molecular_weight * dmole_sources_dX_s[k][j];
      }
      num1 = -cp_dot_species * (h_dot_J * V + h_dot_mole_sources * dV_dxj);
      num2 = (-h_dot_mole_sources * V + rxn.ProblemInfo.Q) * rxn.Thermo->cp(temp_cache, j);
      IJth(J,n_eq, j+1) = (num1 - num2) / cp_dot_species_2;
  }

  // Energy RHS derivative wrt temperature
  num1 = -cp_dot_species * (h_dot_J_T + cp_dot_mole_sources * V);
  num2 = (-h_dot_mole_sources * V + rxn.ProblemInfo.Q) * dcp_dT_dot_species;
  IJth(J,n_eq,n_eq) = (num1 - num2) / cp_dot_species_2;

  return(0);
}

/*
 *:::::::::::::::::::::::::::::::
 * Functions for NLOPT
 *:::::::::::::::::::::::::::::::
 * */

int mass_balance(unsigned int n_species, unsigned int n_atoms, 
                 std::vector<double> Amat, std::vector<double> bvec, 
                 double *x1)
{
    n_atoms += 1;

    double lb[n_species];
    memset(lb, 0.0, sizeof lb);

    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LN_COBYLA, n_species); /* algorithm and dimensionality */
    nlopt_set_lower_bounds(opt, lb);

    my_f_data* fdata = (my_f_data *) malloc(sizeof(my_f_data)); 
    fdata->x1= (double *) malloc(n_species * sizeof(double));

    for (unsigned int k = 0; k < n_species; k++) {
        fdata->x1[k] = x1[k];
    }

    nlopt_set_min_objective(opt, myfunc, fdata);
    my_constraint_data* cdata = (my_constraint_data *) malloc(sizeof(my_constraint_data));

    cdata->b = (double *) malloc(n_atoms * sizeof(double));
    cdata->A = (double *) malloc(n_atoms * n_species * sizeof(double));

    for (unsigned int m = 0; m < n_atoms; m++) {
        cdata->b[m] = bvec[m];
    }

    for (unsigned int kk = 0; kk < n_atoms * n_species; kk++) {
        cdata->A[kk] = Amat[kk];
    }

    double tol[] = {1.0e-08, 1.0e-08, 1.0e-08};
    nlopt_add_equality_mconstraint(opt, n_atoms, myconstraint, cdata, tol);
    nlopt_set_xtol_rel(opt, 1e-4);
    
    for (unsigned int k = 0; k < n_species; k++) {
        if (x1[k] < 0.0) {
           x1[k] = 1.0e-15;  /* some initial guess */
        }
    }

    double minf; /* the minimum objective value, upon return */

    int opt_flag = nlopt_optimize(opt, x1, &minf);

    free(fdata->x1);
    free(fdata);

    free(cdata->b);
    free(cdata->A);
    free(cdata);
    
    nlopt_destroy(opt);

    return(opt_flag);
}

double myfunc(unsigned int n, const double *x, double *grad, void *my_func_data) 
{
    my_f_data *d = (my_f_data *) my_func_data;
    double diff[n];
    double result = 0.0;
    for (unsigned int i = 0; i < n; i++) {
        diff[i] = x[i] - d->x1[i];
        result += diff[i] * diff[i];
    }
    if (grad) {
       for (unsigned int i = 0; i < n; i++) {
           grad[i] = 2.0 * diff[i];
       }
    }
    
    return result;
}

void myconstraint(unsigned int m, double *result, unsigned int n, const double* x, double* grad, void* f_data)
{
    my_constraint_data *d = (my_constraint_data *) f_data;
    if (grad) {
       for (unsigned int i = 0; i < m*n; i++) {
           grad[i] = d->A[i];
       }
    }
    for (unsigned int i = 0; i < m; i++) {
        result[i] = 0.0;
        for (unsigned int j = 0; j < n; j++) {
            result[i] += d->A[i*n + j] * x[j];
        }
        result[i] -= d->b[i];
    }
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

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
