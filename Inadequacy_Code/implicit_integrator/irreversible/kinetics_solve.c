/*
 * -----------------------------------------------------------------
 * $Revision: 4834 $
 * $Date: 2016-12-18 16:07:00 -0700 (Sun, 18 Dec 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): David Sondak @i ICES
 * -----------------------------------------------------------------
 * Description:
 * 
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <fstream>
#include <algorithm>

// user defined functions
#include "reaction_info.h"
#include "inadequacy_model.h"
#include "write_data.h"

/* Header files with a description of contents used */
#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <nlopt.h>                   /* nlopt optimization package */

// Antioch
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
#define RTOL  RCONST(5.0e-05)  /* scalar relative tolerance            */
#define ATOL  RCONST(1.0e-11)  /* scalar relative tolerance            */
#define T0    RCONST(0.0)      /* initial time           */
#define T1    RCONST(0.000001) /* first output time      */
#define DTOUT RCONST(0.000001) /* output time factor     */
#define NOUT  15000            /* number of output times */


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

int kinetics_forward(std::vector<double>& initial_condition, 
             std::vector<double>& time_points, 
             reaction_info* rxn, 
             std::vector<double>& solution)
{
  realtype reltol, t;
  N_Vector y, abstol;
  void *cvode_mem;
  int flag;

  y = abstol = NULL;
  cvode_mem = NULL;

  const unsigned int n_eq = initial_condition.size();

  /* Create serial vector of length n_eq for I.C. and abstol */
  y = N_VNew_Serial(n_eq);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(n_eq); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

  /* Set the scalar relative tolerance */
  reltol = RTOL;

  /* Set the vector absolute tolerance */
  for (unsigned int i = 0; i < n_eq; i++) {
      Ith(abstol,i+1) = ATOL;
  }

  /* Initialize solution field */
  for (unsigned int k = 1; k <= n_eq; k++) {
      Ith(y,k) = initial_condition[k-1];
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

  flag = CVodeSetMaxNumSteps(cvode_mem, 10000000);
  if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

  /* Call CVDense to specify the CVDENSE dense linear solver */
  flag = CVDense(cvode_mem, n_eq);
  if (check_flag(&flag, "CVDense", 1)) return(1);

  /* Set the Jacobian routine to Jac (user-supplied) */
  flag = CVDlsSetDenseJacFn(cvode_mem, Jac);
  if (check_flag(&flag, "CVDlsSetDenseJacFn", 1)) return(1);

  /* Set user data */
  flag = CVodeSetUserData(cvode_mem, rxn);
  if (check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  realtype Tig     = 0.0; // ignition temperature
  realtype time_ig = 0.0; // ignition time

  /* In loop, call CVode, test for error, and test for ignition.
     Break out of loop when ignition conditions have been reached.  */

  //t = 0.0;
  for (unsigned int iout = 1; iout <= NOUT; iout++) {

    realtype tout = iout * DTOUT;

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    if (check_flag(&flag, "CVode", 1)) throw flag;

    //if (Ith(y,n_eq) >= 1500.0) {
    if (Ith(y,1) <= 0.99*initial_condition[0]) {
       time_ig = t;
       Tig = Ith(y,n_eq);
       break;
    }

  }

/*=========================================
 *
 * Determine time shift
 *
 *=========================================*/
  double delta_tig = rxn->ProblemInfo.time_ig - time_ig;
  int n_samp = rxn->ProblemInfo.n_times;

  /* Reinitialize solution field */
  for (unsigned int k = 1; k <= n_eq; k++) {
      Ith(y,k) = initial_condition[k-1];
  }

  //t = 0.0;  // Reset time

  flag = CVodeReInit(cvode_mem, T0, y);
  if (check_flag(&flag, "CVodeReInit", 1)) return(1);

  /* In loop, call CVode and test for error.
     Break out of loop when all samples have been taken.  */
  for (unsigned int iout = 0; iout < n_samp; iout++) {

    realtype tout = time_points[iout] - delta_tig;

    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    if (check_flag(&flag, "CVode", 1)) throw flag;

    /* Store solution */
    for (unsigned int j = 0; j < n_eq; j++) {
        solution[iout * n_eq + j] = Ith(y,j+1);
    }

  }

  /* Ignition data */
  solution[n_eq * n_samp] = Tig;

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
  const unsigned int n_inad        = rxn.InadInfo.n_inad;
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
  // Linearly interpolate the temperature to get value at current t
  double temperature = rxn.LinInterp.Interpolate(t);
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
  const unsigned int n_inad        = rxn.InadInfo.n_inad;
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
  // Linearly interpolate the temperature to get value at current t
  double temperature = rxn.LinInterp.Interpolate(t);
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
  }

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
