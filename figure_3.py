import matplotlib.pyplot as plt
from scipy.signal import argrelmax

from solver.solver import *

parameters = default_parameters()
pulse_amplitude = [-0.6, -1.3, -1.9, -2.8]

t, state = solve(parameters, condition=0, tmax=2000)
pks = argrelmax(state[:, 0])[0]

pk0 = pks[2]
base_cycle = t[pks[3]] - t[pks[2]]
stim_time_list = np.linspace(t[pks[2]], t[pks[3]], 50)

for amp in pulse_amplitude:
    print(amp)
    period = []
    for stim_time in stim_time_list:
        all_t = []
        all_y = []
        times = [stim_time, 50, 2000 - stim_time - 50]
        ic = [-1.26312652e+01, 6.73685294e-01, 5.28447665e-01, 9.60815694e-01, 4.83891944e-07, 3.70080101e-02,
              1.24208141e-01]
        for ix, t_length in enumerate(times):
            if ix == 1:
                parameters[-1] = amp
            else:
                parameters[-1] = 0
            [partial_t, partial_y] = solve(parameters, tmax=t_length, state0=ic)
            all_t.append(partial_t)
            all_y.append(partial_y[:, 0])

            ic = partial_y[-1, :]

        all_t = np.concatenate((all_t[0], all_t[1] + all_t[0][-1], all_t[2] + all_t[0][-1] + all_t[1][-1]))
        all_y = np.concatenate((all_y[0], all_y[1], all_y[2]))

        pks = argrelmax(all_y)[0]
        period.append(all_t[pks[3]] - all_t[pks[2]])

    plt.plot(stim_time_list, 100 * (period - base_cycle) / base_cycle)
