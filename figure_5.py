import matplotlib.pyplot as plt

from solver.solver import *

parameters = default_parameters('dep')
IC = [-4.15794496e+01, 2.51895863e-01, 6.80873415e-02, 6.40403785e-01,
      1.55212288e-04, 1.71282992e-01, 2.76027028e-02]
t, state = solve(parameters, condition=0, tmax=2000, state0=IC)
pks = argrelmax(state[:, 0])[0]

base_cycle = t[pks[1]] - t[pks[0]]
thresh = state[:, 0].mean()

cycle_fraction = [0.4, 0.44, 0.48]
crossings = (state[1:, 0] >= thresh) * (state[:-1, 0] < thresh)
crosses = np.where(crossings)[0]

#cycle_start = t[crosses[0]]
cycle_start = t[pks[0]]-40
for ix, cf in enumerate(cycle_fraction):
    stim_time = cycle_start + cf * base_cycle

    print(stim_time)
    amp = -2
    width = 50

    t, y = perturb(stim_time, parameters, amp, width, ic=IC, tmax=3000)
    plt.subplot(311 + ix)
    plt.plot(t, y)

    plt.plot([stim_time, stim_time], [-50, 0], 'k--')
