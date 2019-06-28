import matplotlib.pyplot as plt

from solver.solver import *

plt.figure()
vclamp = np.arange(-60, -10, 0.1)
currentarray = np.zeros(vclamp.shape)
parameters = default_parameters()
for ix, v in enumerate(vclamp):
    state0 = [v, 0, 0, 0, 0, 0, 0]
    t, state2 = solve(parameters, condition=0, tmax=10000, dt=0.2, clamp_ix=0, state0=state0)
    currentarray[ix] = total_current(state2[-1, :], 0, parameters)

plt.plot(vclamp, currentarray)
plt.xlabel('Voltage (mV)')
plt.ylabel('Current (uA/cm^2)')
plt.title('Control')
