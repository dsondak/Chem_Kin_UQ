#!/bin/py

import numpy as np
import matplotlib.pyplot as plt
import h5py
import config

np.set_printoptions(precision=4)

plt.rcParams.update(config.pars)

fsol_detailed = "detailed_solution.h5"
fsol_reduced = "reduced_solution.h5"
fsol_inad = "inad_solution.h5"

sol_detailed = h5py.File(fsol_detailed, 'r')
sol_reduced = h5py.File(fsol_reduced, 'r')
sol_inad = h5py.File(fsol_inad, 'r')

time_d = sol_detailed['Scenario1/time']
time_r = sol_reduced['Scenario1/time']
time_i = sol_inad['Scenario1/time']

data_d = sol_detailed['Scenario1/truth']
data_r = sol_reduced['Scenario1/truth']
data_i = sol_inad['Scenario1/truth']

t = np.array(time_d)
#t_samps = t[np.where((t >= 0.02104) & (t <= 0.02130))]
#print(t[21001])
#t_samps = t_samps[::10]

#t_samps = np.array([0.02105, 0.021091, 0.0211, 0.021102, 0.021105, 0.021106, 
#                    0.021107,  0.021108, 0.021109, 0.02111, 0.0212])
t_samps = np.array([0.0209, 0.021001, 0.02105, 0.021091, 0.0211, 0.021102, 0.021105, 0.021106, 
                    0.021107,  0.021108, 0.021109, 0.02111, 0.02113, 0.0212, 0.0214])

#[0.02105 0.02109 0.0211 0.021102 0.021104 0.021106 0.021107 0.021108 0.021109 0.02111  0.0212]
#[0.02105          0.0211 0.021102          0.021106 0.021107 0.021108 0.021109 0.02111  0.0212  ]

idxs = np.isin(t, t_samps)
idxs = np.where(idxs)
idxs = idxs[0]
#print(idxs)
#print(t[idxs])
#print(np.size(idxs), np.size(t_samps))

#print(idxs)
#print(data_d[idxs, 0])

Ns = 7

fig_alld, ax_alld = plt.subplots(1,1)
fig_allr, ax_allr = plt.subplots(1,1)
fig_alli, ax_alli = plt.subplots(1,1)
for i in range(Ns):

    fig_s, ax_s = plt.subplots(1,1)
    
    ax_s.plot(time_d, data_d[:,i], ls='-', label=r'Detailed')
    ax_s.plot(time_r, data_r[:,i], ls='--', label=r'Reduced')
    ax_s.plot(time_i, data_i[:,i], ls='-.', label=r'Inadequacy')

    ax_s.plot(t_samps, data_d[idxs,i], ls=' ', marker='^', ms=12, mfc='w', mec='k', mew=3)

    ax_alld.plot(time_d, data_d[:,i], ls='-')
    ax_allr.plot(time_r, data_r[:,i], ls='--')
    ax_alli.plot(time_i, data_i[:,i], ls='-.')
    
    ax_s.legend(loc='best')


fig_T, ax_T = plt.subplots(1,1)
T = data_d[:,-1]
ax_T.plot(time_d, T)

fig_v, ax_v = plt.subplots(1,1)

ax_v.plot(time_i, data_i[:,7], ls='-.', label=r'$H^{\prime}$')
ax_v.plot(time_i, data_i[:,8], ls='-.', label=r'$O^{\prime}$')

ax_v.legend(loc='best')


plt.show()
