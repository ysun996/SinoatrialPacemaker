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


def currentODE(state, t, parameters):
    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    ENa = parameters[5]
    ESi = parameters[6]
    EK = parameters[7]
    #Cm = parameters[8]
    I = parameters[9]

    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]

    Isi = 0.78 * (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m ** 3) * h + gnNa) * (v - ENa)
    Ix1 = 0.8 * x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = ((y ** 2) * gfNa) * (v - ENa) + (v - EK) * (y ** 2) * (-120 / EK)

    Ii = Isi + Ina + Ix1 + Ik1 + If + I

    retcurrent = [Isi, Ina, Ix1, Ik1, If, Ii]

    return retcurrent

#####Pacemaker current at 52%

period = pointdiff[3]
currents = [[0, 0, 0, 0, 0, 0]]
stimtime = period + per[1] * 0.29
stimmax = stimtime + 50

t1 = np.arange(0, stimtime, 0.2)
t2 = np.arange(stimtime, stimmax, 0.2)
t3 = np.arange(stimmax, 5000, 0.2)

parameters1 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]
parameters2 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, -2]
parameters3 = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]

state1 = odeint(PacemakerODE, state0, t1, args=(parameters1,), hmax=0.2)
state2 = odeint(PacemakerODE, state1[-1, :], t2, args=(parameters2,), hmax=0.2)
state3 = odeint(PacemakerODE, state2[-1, :], t3, args=(parameters3,), hmax=0.2)

allstate = np.concatenate([state1, state2, state3])
alltime = np.concatenate([t1, t2, t3])

for states in allstate:
    currentval = [currentODE(states, t, parameters)]
    currents = np.concatenate([currents, currentval], axis = 0)
plt.legend(['Isi', 'Ina', 'Ix1', 'Ik1', 'If'])
for i in range(len(currents)):
    plt.plot(alltime, currents[1:, i])