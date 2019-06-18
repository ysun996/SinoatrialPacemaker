from commonDep import *

phasedata = np.arange(0.3, 0.5, 0.01)
stimcurrent = np.arange(-3, -2, 0.01)
for i in phasedata:
    annipulsefun, timepulse = pulsefun(state0, periodhyp, perhyp, i, 0.7, -2, condition)
    plt.plot(timepulse, annipulsefun[:, 0], label = str(i))
plt.title('Hyperpolarized model - depolarizing current')
plt.ylabel('Voltage (mV)')
plt.legend()