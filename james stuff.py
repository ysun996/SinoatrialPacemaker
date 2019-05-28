import numpy as np
import matplotlib.pyplot as plt
import bisect
from scipy.integrate import odeint


#####Sodium current gating constants
def dzdt(a, b, z):
    hill = a / (a + b)
    timeC = 1 / (a + b)
    rVal = (hill - z) / timeC
    return rVal


def am(Vm):
    max_a_m = 1
    a_m = max_a_m * ((Vm + 47) / (1 - np.exp(-(Vm + 47) / 10)))

    return a_m


def an(Vm):
    a_n = 0.01 * (Vm + 55) / (1 - np.exp(-(Vm + 55 / 10)))
    return a_n


def As(v):
    a_s = 0.05 * (v + 28) / (1 - np.exp(-(v + 28) / 5))
    return a_s


def bs(v):
    b_s = 0.00025 * np.exp(-(v + 28) / 14.93)
    return b_s


def bn(Vm):
    bn = 0.125 * np.exp((-Vm + 65) / 80)
    return bn


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
    max_a_y = 0.00005
    a_y = (max_a_y * np.exp(-(Vm + 32) / 14.93)) / 3
    # a_y = np.exp(-0.0220741*(Vm + 386.9))
    return a_y


def by(Vm):
    max_b_y = 0.001
    b_y = (max_b_y * (Vm + 32) / (1 - np.exp(-(Vm + 32) / 5))) / 3
    # b_y = np.exp(0.052*(Vm - 73.08))
    return b_y


#####Time-independent background potassium current (I_k1)

def i_k1(Vm):
    IPot = 1.3 * (np.exp((Vm + 110) / 25) - 1) / (np.exp((Vm + 60) / 12.5) + np.exp((Vm + 60) / 25)) + \
           0.26 * (Vm + 30) / (1 - np.exp(-(Vm + 30) / 25))
    #    if isinstance(Vm,float)==True:
    #        if Vm > -30.5 and Vm < -29.5:
    #            IPot=1.3*(np.exp((Vm+110)/25)-1)/(np.exp((Vm+60)/12.5)+np.exp((Vm+60)/25))
    #        else:
    #            IPot=1.3*(np.exp((Vm+110)/25)-1)/(np.exp((Vm+60)/12.5)+np.exp((Vm+60)/25))+\
    #                 0.26*(Vm+30)/(1-np.exp(-(Vm+30)/25))
    #    else:
    #        Length=len(Vm)
    #        IPot=np.zeros(Length)
    #        for i in range(Length):
    #            if Vm[i] > -30.5 and Vm[i] < -29.5:
    #                IPot[i]=1.3*(np.exp((Vm[i]+110)/25)-1)/(np.exp((Vm[i]+60)/12.5)+np.exp((Vm[i]+60)/25))
    #            else:
    #                IPot[i]=1.3*(np.exp((Vm[i]+110)/25)-1)/(np.exp((Vm[i]+60)/12.5)+np.exp((Vm[i]+60)/25))+0.26*(Vm[i]+30)/(1-np.exp(-(Vm[i]+30)/25))
    return IPot


# Parameters
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

# Time
tmax = 1000
dt = 0.4
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
    I = parameters[9]
    #    for x in range(1,7):
    #        if state[x] > 1:
    #            state[x]=1

    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]

    #    zArray=[0,d,f,m,h,x1,y]
    Isi = (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m ** 3) * h + gnNa) * (v - ENa)
    Ix1 = x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    # Is = s*(0.6*np.exp((v+110)/25)-1)/(np.exp((v+60)/12.5)+np.exp((v+60)/25))
    If = (y ** 2 * gfNa) * (v - ENa) + (v - EK) * y ** 2 * (-120 / EK)

    Ii = Isi + Ina + Ix1 + Ik1 + If + I
    # print(Ix1)
    dv = -Ii / Cm

    #    rstate = [dv, 0, 0, 0, 0, 0, 0]
    #    aVals=[0,ad(v),af(v),am(v),ah(v),ax1(v),ay(v)]
    #    bVals=[0,bd(v),bf(v),bm(v),bh(v),bx1(v),by(v)]
    #    for x in range(1,7):
    #        rstate[x]=dzdt(aVals[x],bVals[x],zArray[x])
    dd = (ad(v) * (1 - d)) - (bd(v) * d)
    df = (af(v) * (1 - f)) - (bf(v) * f)
    dm = (am(v) * (1 - m)) - (bm(v) * m)
    dh = (ah(v) * (1 - h)) - (bh(v) * h)
    dx1 = (ax1(v) * (1 - x1)) - (bx1(v) * x1)
    # ds = As(v) * (1-s) - bs(v) * s
    dy = (ay(v) * (1 - y)) - (by(v) * y)
    rstate = [dv, dd, df, dm, dh, dx1, dy]

    return rstate


def currentVals(state, t, parameters):
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
    Ina = (gNa * (m ** 3) * h + gnNa) * (v - ENa)
    Ix1 = x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = y ** 2 * (gfNa * (v - ENa) + (v - EK) * (-120 / EK))
    # If = y**2*(4.6 * (v - ENa) + (v-EK) * 7.4)
    # Is = s*(0.6*np.exp((v+110)/25)-1)/(np.exp((v+60)/12.5)+np.exp((v+60)/25))
    Ii = Isi + Ina + Ix1 + Ik1 + If + I
    return [Isi, Ina, Ix1, Ik1, If, Ii]


state0 = [-50, 0, 0, 0, 0, 0, 0]

pulseYN = False
if pulseYN:
    for x in range(20):
        if x == 0:
            tIndexL = bisect.bisect_left(t, 20)
            t1 = t[0:tIndexL]
            tplot = t1
            state = odeint(PacemakerODE, state0, t1, args=(parameters,))
            stateFinal = state
        elif x > 0:
            CycleValue = 20 * (x + 1)
            stimCycleValue = 0.05 + 20 * x
            tIndexR = bisect.bisect_left(t, CycleValue)
            tIndexStim = bisect.bisect_left(t, stimCycleValue)
            parameters[9] = 1000
            t1 = t[tIndexL:tIndexStim]
            tplot = np.append(tplot, t1)
            state1 = state[-1, :]
            state = odeint(PacemakerODE, state1, t1, args=(parameters,))
            stateFinal = np.concatenate((stateFinal, state), axis=0)
            parameters[9] = 0
            t1 = t[tIndexStim:tIndexR]
            tplot = np.append(tplot, t1)
            state1 = state[-1, :]
            state = odeint(PacemakerODE, state1, t1, args=(parameters,))
            stateFinal = np.concatenate((stateFinal, state), axis=0)
            tIndexL = tIndexR

parameters[9] = 0
tmax = 2000
dt = 0.4
t = np.arange(0, tmax, dt)
#############################################################################################
# v vs. dv
state = odeint(PacemakerODE, state0, t, args=(parameters,))
plt.plot(t, state[:,0])
# dstate = [[0, 0, 0, 0, 0, 0, 0]]
# for w in range(len(state[:, 0])):
#     d1state = [PacemakerODE(state[w, :], t, parameters)]
#     dstate = np.concatenate((dstate, d1state), axis=0)
# plt.plot(state[:, 0], dstate[1:, 0])
#############################################################################################
# parameters[9] = 0
# state = odeint(PacemakerODE, state0, t, args=(parameters,))
# dstate=[[0,0,0,0,0,0]]
# for w in range(len(state[:,0])):
#    d1state=[currentVals(state[w,:],t,parameters)]
#    dstate=np.concatenate((dstate,d1state),axis=0)


# Ilabels=['Isi','Ina','Ix1','Ik1','If','Ii']
# for w in range(0,5):
#    plt.plot(t,dstate[1:,w],label=Ilabels[w])
# plt.legend()
# zlabels=['v','d','f','m','h','x1','y']
# for w in range(1,7):
#    plt.plot(state[:,0],state[:,w],label=zlabels[w])
# plt.plot(t,state[:,w],label=zlabels[w])
# voltage = np.arange(-70, 30, 0.2)

Ibias = np.arange(0, 10, 0.5)

# plt.plot(t,state[:,0])
# for v in range(len(Ibias)):
#    parameters[9] = Ibias[v]
#    state = odeint(PacemakerODE, state0, t, args=(parameters,))
#    plt.subplot(5,5,v+1)
#    plt.plot(t,state[:,0])

plt.ylabel('dV/dt')
plt.xlabel('V')