import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft
from scipy.signal import argrelextrema
from common import *
from commonDep import *


state0 = [0, 0, 0, 0, 0, 0, 0]

state = odeint(PacemakerODE, state0, t, args=(parameters,))
pointdiff = argrelextrema(np.array(state[:,0]), np.greater)[0] * 0.2
per = np.diff(pointdiff)


#######Annihilation stuff
period = pointdiff[3]
phasedata = np.arange(0.29, 0.3, 0.01)
stimcurrent = np.arange(-3.5, -2.5, 0.01)
# plt.subplot(3, 1, 1)


#Parameter sweep for current
for i in stimcurrent:
    stimtime = period + per[1] * 0.29
    stimmax = stimtime + 50

    t1 = np.arange(0, stimtime, 0.2)
    t2 = np.arange(stimtime, stimmax, 0.2)
    t3 = np.arange(stimmax, 5000, 0.2)

    ###Parametersa
    # for j in stimcurrent:
    parameters1 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]
    parameters2 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, i]
    parameters3 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]

    state1 = odeint(PacemakerODE, state0, t1, args=(parameters1,), hmax=0.2)
    state2 = odeint(PacemakerODE, state1[-1, :], t2, args=(parameters2,), hmax=0.2)
    state3 = odeint(PacemakerODE, state2[-1, :], t3, args=(parameters3,), hmax=0.2)

    truestate = np.concatenate([state1[:, 0], state2[:, 0], state3[:,0]])
    truetime = np.concatenate([t1, t2, t3])
    plt.plot(truetime, truestate, label = str(i))
    # plt.figure()
plt.legend()
plt.title('Depressed model')
plt.ylabel('Voltage (mV)')
plt.xlabel('Time (ms)')