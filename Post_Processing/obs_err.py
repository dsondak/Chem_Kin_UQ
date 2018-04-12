#!/bin/py

import numpy as np
import matplotlib.pyplot as plt
import h5py
import config

plt.rcParams.update(config.pars)

# Read in solutions
fsol_detailed = "detailed_solution.h5"

sol_detailed = h5py.File(fsol_detailed, 'r')

time_d = sol_detailed['Scenario1/time']

data_d = sol_detailed['Scenario1/truth']

# Now get some samples from the truth data
t = np.array(time_d)
t_samps = np.array([0.0209, 0.021001, 0.02105, 0.021091, 0.0211, 0.021102, 0.021105, 0.021106, 
                    0.021107,  0.021108, 0.021109, 0.02111, 0.02113, 0.0212, 0.0214])
idxs = np.isin(t, t_samps)
idxs = np.where(idxs)
idxs = idxs[0]

obs_data = np.zeros([np.size(t_samps), 10])
#for ii in range(10):
#    obs_data[:,ii] = data_d[idxs,ii]
for ii in range(10):
    for jj in range(np.size(t_samps)):
        mean_data = data_d[idxs[jj],ii]
        pert = np.sqrt(0.0025 * mean_data * mean_data) * np.random.randn()
        obs_data[jj,ii] = mean_data + pert

ig_data = sol_detailed['Scenario1/ignition']
#print(ig_data[0], ig_data[1])

# Create truth_data.h5
td = h5py.File('truth_test.h5', 'w')
grp = td.create_group("Scenario1")
grp.create_dataset("ignition", data=ig_data)
grp.create_dataset("time", data=t_samps)
dset = grp.create_dataset("truth", data=obs_data)
dset.attrs["Scenario Parameters"] = np.array([1.0, 5.0e+06])

Ns = 7

fig_alld, ax_alld = plt.subplots(1,1)
for i in range(Ns):

    fig_s, ax_s = plt.subplots(1,1)
    
    ax_s.plot(time_d, data_d[:,i], ls='-', label=r'Detailed')
    ax_s.plot(t_samps, obs_data[:,i], ls=' ', marker='^', ms=10, mfc='w',
mec='k', mew=1.5)

    conc_mean = data_d[:,i]
    
    conc_rand = np.sqrt(0.01 * conc_mean * conc_mean) * np.random.randn()
    conc_rand_up = conc_mean + conc_rand
    conc_rand_lo = conc_mean - conc_rand

    ax_s.plot(time_d, conc_rand_up, ls='-', lw=3, color='#D16103')
    ax_s.plot(time_d, conc_rand_lo, ls='-', lw=3, color='#D16103')

    ax_s.fill_between(time_d, conc_rand_up, conc_rand_lo, color='#D16103',
alpha=0.35)

    ax_alld.plot(time_d, data_d[:,i], ls='-')
    
    ax_s.legend(loc='best')

ax_alld.plot(np.array([ig_data[0], ig_data[0]]), np.array([0.0, 2.0]), ls='-',
lw=2)


fig_T, ax_T = plt.subplots(1,1)
T = data_d[:,-1]
ax_T.plot(time_d, T)
ax_T.plot(t_samps, obs_data[:,-1], ls=' ', marker='^', ms=12, mfc='w', mec='k',
mew=1.5)

plt.show()
