import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft
from scipy.signal import argrelextrema
#comment

#####Sodium current gating constants
def am(Vm):
    max_a_m = 1

    a_m = max_a_m * ((Vm + 47) / (1 - np.exp(-(Vm + 47) / 10)))
    return a_m

def bm(Vm):
    max_b_m = 40
    b_m = max_b_m * np.exp(-(Vm + 72) / 17.86)


    return b_m

def ah(Vm):
    max_a_h = 0.0085
    a_h = max_a_h * np.exp(-(Vm + 71) / 5.43)


    return a_h

def bh(Vm):
    max_b_h = 2.5
    b_h = max_b_h / (1 + np.exp(-(Vm + 10) / 12.2))

    return b_h

#####Slow inward gating constants

def ad(Vm):
    max_a_d = 0.005
    a_d = max_a_d * ((Vm + 34) / (1 - np.exp(-(Vm + 34) / 10)))
    return a_d

def bd(Vm):
    max_b_d = 0.05
    b_d = max_b_d * (np.exp(-(Vm + 34) / 6.67))

    return b_d

def not_d(Vm):
    ### d' variable needed for the slow inward current
    not_d = 1 / (1 + np.exp(-(Vm + 15) / 6.67))

    return not_d

def af(Vm):
    max_a_f = 0.002468
    a_f = max_a_f * np.exp(-(Vm + 47) / 20)

    return a_f

def bf(Vm):
    max_b_f = 0.05
    b_f = max_b_f / (1 + np.exp(-(Vm + 13) / 11.49))


    return b_f

##### Constants for delayed rectifier (I_x1)

def ax1(Vm):
    max_a_x1 = 0.0025
    a_x1 = (max_a_x1 * np.exp((Vm + 30) / 12.11)) / (1 + np.exp((Vm + 30) / 50))


    return a_x1

def bx1(Vm):
    max_b_x1 = 0.0065
    b_x1 = max_b_x1 * (np.exp(-(Vm - 20) / 16.67)) / (1 + np.exp(-(Vm - 20) / 25))

    return b_x1

def l_bar_x1(Vm):
    l_bar = 2.25 * ((np.exp((Vm + 95) / 25) - 1) / (np.exp((Vm + 45) / 25)))

    return l_bar

#####Constants for pacemaker current (I_y)

def ay(Vm):
    max_a_y = 0.00005/3.0811907104922107
    a_y = max_a_y * np.exp(-(Vm + 14.25) / 14.93)

    return a_y

def by(Vm):
    max_b_y = 0.001/3.0811907104922107
    b_y = max_b_y * (Vm + 14.25) / (1 - np.exp(-(Vm + 14.25) / 14.93))

    return b_y


#####Time-independent background potassium current (I_k1)

def i_k1(Vm):
    IPot = 1.3 * (np.exp((Vm + 110) / 25) - 1) / (np.exp((Vm + 60) / 12.5) + np.exp((Vm + 60) / 25)) + \
           0.26 * (Vm + 30) / (1 - np.exp(-(Vm + 30) / 25))
    return IPot



#Parameters
gNa = 4.4
gnNa = 0.066
gSi = 0.5175
gnSi = 0.161
gfNa = 1.2
ENa = 40
ESi = 70
EK = -93
Cm = 6
Ibias = -2
Istim = 0
parameters = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, Ibias, Istim]

#Time
tmax = 10000
dt = 0.2
t = np.arange(0, tmax, dt)

def PacemakerODE(state, t, parameters):
    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    ENa = parameters[5]
    ESi = parameters[6]
    EK = parameters[7]
    Cm = parameters[8]
    Ibias = parameters[9]
    Istim = parameters[10]

    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]

    Isi = 0.78 * (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m**3) * h + gnNa) * (v - ENa)
    Ix1 = 0.8 * x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = ((y**2) * gfNa) * (v - ENa) + (v-EK) * (y**2) * (-120/EK)

    Ii = Isi + Ina + Ix1 + Ik1 + If + Ibias + Istim

    dv = -Ii/Cm

    dd = (ad(v) * (1-d)) - (bd(v) * d)
    df = (af(v) * (1-f)) - (bf(v) * f)

    dm = (am(v) * (1-m)) - (bm(v) * m)
    dh = (ah(v) * (1-h)) - (bh(v) * h)

    dx1 = (ax1(v) * (1-x1)) - (bx1(v) * x1)
    #print(dx1)

    dy = (ay(v) * (1-y)) - (by(v) * y)

    rstate = [dv, dd, df, dm, dh, dx1, dy]

    return rstate

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


######Currents

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

######Phase-plane relationship ######

dv = []

for states in state1:
    pacemaker1 = PacemakerODE(states, t1, parameters1)
    dv.append(pacemaker1[0])


for states in state2:
    pacemaker2 = PacemakerODE(states, t2, parameters2)
    dv.append(pacemaker2[0])



for states in state3:
    pacemaker3 = PacemakerODE(states, t3, parameters3)
    dv.append(pacemaker3[0])


# allstate = np.concatenate([pacemaker1, pacemaker2, pacemaker3])
plt.plot(allstate[:, 0], dv)
plt.title('Phase-plane relationship at 29%')






