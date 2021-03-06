from updated_files.common import *
from commonHyp import *


state0 = [0, 0, 0, 0, 0, 0, 0]
state = odeint(PacemakerODE, state0, t, args=(parameters,))
pointdiff = argrelextrema(np.array(state[:,0]), np.greater)[0] * 0.2
per = np.diff(pointdiff)

#######Annihilation stuff
period1 = pointdiff[3]
phasedata = np.arange(0.3, 0.5, 0.01)
goodphases = []
stimcurrent = np.arange(-3, -2, 0.01)
plt.subplot(2,1,1)
for i in phasedata:
    stimtime = period1 + per[1] * i
    stimmax = stimtime + 50

    t1 = np.arange(0, stimtime, 0.2)
    t2 = np.arange(stimtime, stimmax, 0.2)
    t3 = np.arange(stimmax, 5000, 0.2)

    ###Parametersa
    parameters1 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]
    parameters2 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, -2]
    parameters3 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]

    state1 = odeint(PacemakerODE, state0, t1, args=(parameters1,))
    state2 = odeint(PacemakerODE, state1[-1, :], t2, args=(parameters2,))
    state3 = odeint(PacemakerODE, state2[-1, :], t3, args=(parameters3,))

    truestate = np.concatenate([state1[:, 0], state2[:, 0], state3[:,0]])
    truem = np.concatenate([state1[:,1], state2[:, 1], state3[:,1]])
    truetime = np.concatenate([t1, t2, t3])

    timediff = argrelextrema(truestate, np.greater)
    timeper = timediff[0] * 0.2
    difference = np.diff(timeper)
    plt.plot(truetime, truestate, label = str(int((i*100))) + '%')
plt.figure()
plt.gca().axes.get_xaxis().set_visible(False)
plt.title('Hyperpolarized model - depolarizing current')
plt.ylabel('Voltage (mV)')
plt.legend()

#######Depolarizing current
period2 = pointdiff[3]
phasedata = np.arange(0.75, 0.95, 0.01)
stimcurrent = np.arange(-3, -2, 0.01)
plt.subplot(2,1,2)
for i in phasedata:
    stimtime = period2 + per[1] * i
    stimmax = stimtime + 50

    t1 = np.arange(0, stimtime, 0.2)
    t2 = np.arange(stimtime, stimmax, 0.2)
    t3 = np.arange(stimmax, 5000, 0.2)

    ###Parametersa
    # for j in stimcurrent:
    parameters1 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]
    parameters2 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, 2]
    parameters3 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]

    state1 = odeint(PacemakerODE, state0, t1, args=(parameters1,))
    state2 = odeint(PacemakerODE, state1[-1, :], t2, args=(parameters2,))
    state3 = odeint(PacemakerODE, state2[-1, :], t3, args=(parameters3,))

    truestate = np.concatenate([state1[:, 0], state2[:, 0], state3[:,0]])
    truem = np.concatenate([state1[:,1], state2[:, 1], state3[:,1]])
    truetime = np.concatenate([t1, t2, t3])

    timediff = argrelextrema(truestate, np.greater)
    timeper = timediff[0] * 0.2
    difference = np.diff(timeper)
    plt.plot(truetime, truestate, label=str(int((i * 100))) + '%')

plt.title('Hyperpolarized model - hyperpolarizing current')
plt.xlabel('Time (ms)')
plt.legend()
plt.figure()










