import matplotlib.pyplot as plt

from solver.solver import *

plt.figure()
vclamp = np.arange(-60, -10, 0.5) + 1e-6  # singularity in some voltage eq'ns
currentarray = np.zeros(vclamp.shape)

for ix, model_type in enumerate(['base', 'dep', 'hyp']):
    state0 = None
    state = None

    plt.subplot(311 + ix)
    parameters = default_parameters(state=model_type)
    for iy, v in enumerate(vclamp):
        if state0 is None:
            state0 = [v, 0.004, 0.85, 0.192, 0.0267, 0.00437, 0.204]
            tmax = 5000
        else:
            state0 = state[-1, :]
            state0[0] = v
            tmax = 300

        t, state = solve(parameters, condition=0, tmax=tmax, dt=0.5, clamp_ix=0, state0=state0)
        currentarray[iy] = total_current(state[-1, :], parameters)

    plt.plot(vclamp, currentarray)
    plt.plot(vclamp, [0] * len(vclamp))
