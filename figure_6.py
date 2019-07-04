import matplotlib.pyplot as plt

from solver.solver import *

plt.figure()
vclamp = np.arange(-60, -10, 0.1)
currentarray = np.zeros(vclamp.shape)

for ix, model_type in enumerate(['base', 'dep', 'hyp']):
    plt.subplot(311 + ix)
    parameters = default_parameters(state=model_type)
    for iy, v in enumerate(vclamp):
        state0 = [v, 0, 0, 0, 0, 0, 0]
        t, state = solve(parameters, condition=0, tmax=5000, dt=0.2, clamp_ix=0, state0=state0)
        currentarray[iy] = total_current(state[-1, :], 0, parameters)

    plt.plot(vclamp, currentarray)
    plt.plot(vclamp, [0] * len(vclamp))
