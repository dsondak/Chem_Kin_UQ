#!/bin/py

import numpy as np
import matplotlib.pyplot as plt
import h5py
import config

plt.rcParams.update(config.pars)

# Here are two species-to-atoms transition matrices
# First detailed
#                 H2 O2 H  O  OH HO2 H2O H2O2
A_d = np.array([ [2, 0, 1, 0, 1, 1,  2,  2], 
                 [0, 2, 0, 1, 1, 2,  1,  2] ])
# Next reduced
#                 H2 O2 H  O  OH HO2 H2O
A_r = np.array([ [2, 0, 1, 0, 1, 1,  2], 
                 [0, 2, 0, 1, 1, 2,  1] ])

# Files with useful data in them
fsol_detailed = "detailed_solution.h5"
fsol_reduced = "reduced_solution.h5"

fsol_detailed_new = "detailed_solution_new.h5"
fsol_reduced_new = "reduced_solution_new.h5"

frhs_detailed = "detailed_rhs.h5"
frhs_reduced = "reduced_rhs.h5"

# Open up those files
sol_detailed = h5py.File(fsol_detailed, 'r')
sol_reduced  = h5py.File(fsol_reduced, 'r')

sol_detailed_new = h5py.File(fsol_detailed_new, 'r')
sol_reduced_new  = h5py.File(fsol_reduced_new, 'r')

rhs_detailed = h5py.File(frhs_detailed, 'r')
rhs_reduced  = h5py.File(frhs_reduced, 'r')

# Probe solution files for some desired data
time_d = sol_detailed['Scenario1/time']
time_r = sol_reduced['Scenario1/time']

time_d_new = sol_detailed_new['Scenario1/time']
time_r_new = sol_reduced_new['Scenario1/time']

data_d = sol_detailed['Scenario1/truth']
data_r = sol_reduced['Scenario1/truth']

data_d_new = sol_detailed_new['Scenario1/truth']
data_r_new = sol_reduced_new['Scenario1/truth']

# Probe RHS files for some stuff
fx_detailed = rhs_detailed['rhs/fx']
fx_reduced  = rhs_reduced['rhs/fx']
fT_detailed = rhs_detailed['rhs/fT']
fT_reduced  = rhs_reduced['rhs/fT']

# Here comes a plot of temperature
figT, axT = plt.subplots(1,1, figsize=(12,9))
axT.plot(time_d, data_d[:,-1], 
        label=r'Detailed')
axT.plot(time_d_new, data_d_new[:,-1], ls='--',
        label=r'Detailed')
axT.plot(time_r, data_r[:,-1], 
        label=r'Reduced')
axT.plot(time_r_new, data_r_new[:,-1], ls='--',
        label=r'Reduced')

axT.legend(loc='best')
axT.set_xlabel(r'$t$')
axT.set_ylabel(r'$T$')

# Here comes a plot of a species (here OH)
figs, axs = plt.subplots(1,1, figsize=(12,9))
axs.plot(time_d, data_d[:,4], 
        label=r'Detailed')
axs.plot(time_r, data_r[:,4], 
        label=r'Reduced')

axs.legend(loc='best')
axs.set_xlabel(r'$t$')
axs.set_ylabel(r'$x_{_{OH}}$')

# Here's a plot of RHS functions for H2
figf, axf = plt.subplots(1,1, figsize=(12,9))
axf.plot(time_d[:-1], fx_detailed[:,0], label='Detailed')
axf.plot(time_d[:-1], fx_reduced[:,0], label='Reduced')

axf.legend(loc='best')
axf.set_xlabel(r'$t$')
axf.set_ylabel(r'$f_{x}$')

# Here's a plot of energy RHS functions
figE, axE = plt.subplots(1,1, figsize=(12,9))
axE.plot(time_d[:-1], fT_detailed, label='Detailed')
axE.plot(time_d[:-1], fT_reduced, label='Reduced')

axE.legend(loc='best')
axE.set_xlabel(r'$t$')
axE.set_ylabel(r'$f_{T}$')

plt.show()
