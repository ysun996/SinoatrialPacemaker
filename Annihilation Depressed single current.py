from commonDep import *

#####Pacemaker current at 52%

currents = [[0, 0, 0, 0, 0, 0]]
stim = pulsefun(period, per, 0.29, -2, -2)

for states in stim[0]:
    currentval = [currentODE(states, t, parametersdep)]
    currents = np.concatenate([currents, currentval], axis = 0)
plt.legend(['Isi', 'Ina', 'Ix1', 'Ik1', 'If'])
for i in range(len(currents)):
    plt.plot(stim[1], currents[1:, i])