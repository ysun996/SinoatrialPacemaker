import matplotlib.pyplot as plt

from solver.solver import *

plt.figure()

"""
Vm
"""

plt.subplot(411)

parameters = default_parameters()
[t, state] = solve(parameters)
plt.plot(t, state[:, 0])
plt.xlim([250, 1250])

"""
Im
"""

plt.subplot(412)

im = np.zeros(len(state))
for ix in np.arange(len(im)):
    im[ix] = total_current(state[ix, :], t, parameters)
plt.plot(t, im)
plt.xlim((250, 1250))

"""
Individual currents
"""

plt.subplot(413)

currents = np.array(list(map(lambda x: current(x, t, parameters), state)))
Ilabels = ['Isi', 'Ina', 'Ix1', 'Ik1', 'If', 'Ii']
for w in range(0, 5):
    plt.plot(t, currents[:, w], label=Ilabels[w])

plt.xlim([250, 1250])
plt.legend(Ilabels)

"""
Phase plane
"""
plt.subplot(414)

plt.plot(state[:, 0], np.gradient(state[:, 0], 0.4))
