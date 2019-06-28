from updated_files.commonDep import *

#######Annihilation stuff
phasedata = np.arange(0.29, 0.3, 0.01)
stimcurrent = np.arange(-3.5, -2.5, 0.01)
# plt.subplot(3, 1, 1)


#Parameter sweep for current
for i in stimcurrent:
    annipulsefun, timepulse = pulsefun(state0, period, per, 0.29, -2, i)
    plt.plot(timepulse, annipulsefun[:, 0], label = str(i))
    # plt.figure()a
plt.legend()
plt.title('Depressed model')
plt.ylabel('Voltage (mV)')
plt.xlabel('Time (ms)')