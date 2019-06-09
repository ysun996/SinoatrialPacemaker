from commonDep import *

#####Pacemaker current at 52%

currents = [[0, 0, 0, 0, 0, 0]]
stim, timepulse = pulsefun(state0, period, per, 0.29, -2, -2)

for states in stim:
    currentval = [currentODE(states, t, parametersdep)]
    currents = np.concatenate([currents, currentval], axis = 0)
plt.legend(['Isi', 'Ina', 'Ix1', 'Ik1', 'If'])
for i in range(len(currents)-1):
    plt.plot(timepulse, currents[1:, i])
