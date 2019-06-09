from commonDep import *

phasedata = np.arange(0.75, 0.95, 0.01)
stimcurrent = np.arange(-3, -2, 0.01)

for i in phasedata:
    annipulsefun, timepulse = pulsefun(state0, periodhyp, perhyp, i, 0.7, 2)
    plt.plot(timepulse, annipulsefun[:, 0], label = str(i))
plt.title('Hyperpolarized model - hyperpolarizing current')
plt.xlabel('Time (ms)')
plt.legend()