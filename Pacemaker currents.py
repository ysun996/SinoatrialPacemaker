import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.fftpack import fft
from scipy.signal import argrelextrema

#####Sodium current gating constants
def am(Vm):
    max_a_m = 1
    # checker = 1 - np.exp(-(Vm + 47) / 10)
    # if type(Vm) != np.float64:
    #     modarray = np.zeros(len(Vm))
    #     for i in range(len(Vm)):
    #         checker = 1 - np.exp(-(Vm[i] + 47) / 10)
    #         if checker != 0:
    #             modarray[i] = max_a_m * ((Vm[i] + 47) / (1 - np.exp(-(Vm[i] + 47) / 10)))
    #         else:
    #             modarray[i] = max_a_m * ((Vm[i] + 0.001 + 47) / (1 - np.exp(-(Vm[i] + 0.001 + 47) / 10)))
    #     return modarray
    # else:
    #     if checker != 0:
    #         a_m = max_a_m * ((Vm + 47) / (1 - np.exp(-(Vm + 47) / 10)))
    #         return a_m
    #     return max_a_m * ((47.001 + 47) / (1 - np.exp(-(47.001 + 47) / 10)))
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
    # checker = 1 - np.exp(-(Vm + 34) / 10)
    # if type(Vm) != np.float64:
    #     modarray = np.zeros(len(Vm))
    #     for i in range(len(Vm)):
    #         checker = 1 - np.exp(-(Vm[i] + 34) / 10)
    #         if checker != 0:
    #             modarray[i] = max_a_d * ((Vm[i] + 34) / (1 - np.exp(-(Vm[i] + 34) / 10)))
    #         else:
    #             modarray[i] = max_a_d * ((Vm[i]+ 0.001 + 34) / (1 - np.exp(-(Vm[i] + 0.001 + 34) / 10)))
    #     return modarray
    # else:
    #     if checker != 0:
    #         a_d = max_a_d * ((Vm + 34) / (1 - np.exp(-(Vm + 34) / 10)))
    #         return a_d
    #     return max_a_d * ((Vm + 0.001 + 34) / (1 - np.exp(-(Vm + 0.001 + 34) / 10)))
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
    # checker = 1-np.exp(-(Vm + 17)/14.93)
    # if type(Vm) != np.float64:
    #     modarray = np.zeros(len(Vm))
    #     for i in range(len(Vm)):
    #         checker = 1-np.exp(-(Vm[i] + 17)/14.93)
    #         if checker != 0:
    #             modarray[i] = max_b_y * (Vm[i] + 17)/(1-np.exp(-(Vm[i] + 17)/14.93))
    #         else:
    #             modarray[i] = max_b_y * (Vm[i] + 0.001 + 17)/(1-np.exp(-(Vm[i] + 0.001 + 17)/14.93))
    #     return modarray
    # else:
    #     if checker != 0:
    #         b_y = max_b_y * (Vm + 17)/(1-np.exp(-(Vm + 17)/14.93))
    #         return b_y
    #     return max_b_y * (Vm + 0.001 + 17)/(1-np.exp(-(Vm + 0.001 + 17)/14.93))
    b_y = max_b_y * (Vm + 14.25) / (1 - np.exp(-(Vm + 14.25) / 14.93))
    # b_y = np.exp(0.052*(Vm - 73.08))
    return b_y

# def y_inf(Vm):
#     yinf = ay(Vm)/(ay(Vm)+by(Vm))
#     return yinf
#
# initcon = np.arange(-100, 0, 0.5)
# plt.plot(initcon, ay(initcon))
# plt.plot(initcon, by(initcon))
# plt.plot(initcon, y_inf(initcon))

#####Time-independent background potassium current (I_k1)

def i_k1(Vm):
    # check = 1-np.exp(-(Vm+30)/25)
    # i_k1pt1 = 1.3 * (
    #         (np.exp((Vm + 110) / 25) - 1) / ((np.exp((Vm + 60) / 12.5)) + (np.exp((Vm + 60) / 25))))
    # if type(Vm) != np.float64:
    #     modarray = np.zeros(len(Vm))
    #
    #     for i in range(len(Vm)):
    #         checker = 1-np.exp(-(Vm[i]+30)/25)
    #         i_k1pt1 = 1.3 * (
    #                 (np.exp((Vm[i] + 110) / 25) - 1) / ((np.exp((Vm[i] + 60) / 12.5)) + (np.exp((Vm[i] + 60) / 25))))
    #         if checker != 0:
    #
    #             i_k1pt2 = 0.26 * ((Vm[i] + 30) / (1 - np.exp(-(Vm[i] + 30) / 25)))
    #             modarray[i] = i_k1pt1 + i_k1pt2
    #         else:
    #
    #             modarray[i] = i_k1pt1 + (0.26*((Vm[i]+30.001)/(1-np.exp(-(Vm[i]+30.001)/25))))
    #     return modarray
    # if check != 0:
    #     i_k1pt2 = 0.26*((Vm+30)/(1-np.exp(-(Vm+30)/25)))
    #     return i_k1pt1 + i_k1pt2
    # return i_k1pt1 + 0.26*((Vm+30.001)/(1-np.exp(-(30.001+30)/25)))
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
I = 0
parameters = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, I]

#Time
tmax = 4000
dt = 0.2
t = np.arange(400, tmax, dt)

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
    I = parameters[9]

    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]

    Isi = (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m**3) * h + gnNa) * (v - ENa)
    Ix1 = x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = ((y**2) * gfNa) * (v - ENa) + (v-EK) * (y**2) * (-120/EK)

    Ii = Isi + Ina + Ix1 + Ik1 + If + I

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

state0 = [-66, 0, 0, 0, 0, 0, 0]

state = odeint(PacemakerODE, state0, t, args=(parameters,))
dv = []
for points in state:
    pacemaker = PacemakerODE(points, t, parameters)
    dv.append(pacemaker[0])
# state1 = state[190]
# plt.plot(state[:,0], state[:,0]/t)
plt.plot(t, state[:,0])
# plt.figure()
#plt.plot(state[:,0], dv)

#######Fast Fourier transform
volts = np.zeros(len(state[:,0]))
for ind in range(len(volts)):
    volts[ind] = state[ind, 0]
CL = argrelextrema(volts, np.greater) ####Cycle length array, numbers refer to indices
timpoints = CL[0] * 0.2




# state1 = odeint(PacemakerODE, state02, t, args=(parameters,),atol=1e-6, rtol=1e-6, hmax=0.05)
# plt.plot(t, state1[:,0])
# plt.plot(state[:,0], state[:,4])
# v = np.arange(-100, 50, 1)
#
# def hODE(h, t, parameters):
#
#     v = parameters
#
#     dh = (ah(v) * (1-h)) - (bh(v) * h)
#
#     return dh
#
# param = np.arange(0, 1, 0.01)
# for i in v:
#     parameters = i
#     stateh = odeint(hODE, 0, t, args=(parameters,))
#     pltpoint = np.full(len(stateh), i)
#
#     plt.plot(pltpoint, stateh)
#print(state[:,0])
# voltage = np.arange(-70, 30, 0.02)
# state2 = odeint(PacemakerODE, state1, t, args=(parameters,))
# plt.plot(t, state2[:,0])
# for i in range(1, 5):
#     plt.plot(voltage, state[:, i])

# def currentODE(state, t, parameters):
#     gNa = parameters[0]
#     gnNa = parameters[1]
#     gSi = parameters[2]
#     gnSi = parameters[3]
#     gfNa = parameters[4]
#     ENa = parameters[5]
#     ESi = parameters[6]
#     EK = parameters[7]
#     #Cm = parameters[8]
#     I = parameters[9]
#
#     v = state[0]
#     d = state[1]
#     f = state[2]
#     m = state[3]
#     h = state[4]
#     x1 = state[5]
#     y = state[6]
#
#     Isi = (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
#     Ina = (gNa * (m**3) * h + gnNa) * (v - ENa)
#     Ix1 = x1 * l_bar_x1(v)
#     Ik1 = i_k1(v)
#     If = ((y**2) * gfNa) * (v - ENa) + (v-EK) * (y**2) * (-120/EK)
#
#     Ii = Isi + Ina + Ix1 + Ik1 + If + I
#
#     retcurrent = [Isi, Ina, Ix1, Ik1, If, Ii]
#
#     return retcurrent
#
# currents = [[0, 0, 0, 0, 0, 0]]
# for states in state:
#     currentval = [currentODE(states, t, parameters)]
#     currents = np.concatenate([currents, currentval], axis = 0)
#
# for i in range(len(currents)):
#     plt.plot(t, currents[1:, i])
# plt.legend(['Isi', 'Ina', 'Ix1', 'Ik1', 'If'])






