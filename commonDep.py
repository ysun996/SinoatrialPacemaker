from common import *
from scipy.integrate import odeint
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt

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

    if Ibias < 0:
        Ii = 0.78 * Isi(gSi, gnSi, d, f, v, ESi) + Ina(gNa, m, h, gnNa, v, ENa) + 0.8 * Ix1(x1, v) + i_k1(v) + \
             If(y, gfNa, ENa, EK, v) + Ibias + Istim
    else:
        Ii = Isi(gSi, gnSi, d, f, v, ESi) + Ina(gNa, m, h, gnNa, v, ENa) + Ix1(x1, v) + i_k1(v) + \
             If(y, gfNa, ENa, EK, v) + Ibias + Istim

    dv = -Ii/Cm

    dd = (ad(v) * (1-d)) - (bd(v) * d)
    df = (af(v) * (1-f)) - (bf(v) * f)

    dm = (am(v) * (1-m)) - (bm(v) * m)
    dh = (ah(v) * (1-h)) - (bh(v) * h)

    dx1 = (ax1(v) * (1-x1)) - (bx1(v) * x1)

    dy = (ay(v) * (1-y)) - (by(v) * y)

    rstate = [dv, dd, df, dm, dh, dx1, dy]

    return rstate

def currentODE(state, t, parameters):
    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    ENa = parameters[5]
    ESi = parameters[6]
    EK = parameters[7]
    I = parameters[9]

    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]

    Is = 0.78 * Isi(gSi, gnSi, d, f, v, ESi)
    In = Ina(gNa, m, h, gnNa, v, ENa)
    Ix = 0.8 * Ix1(x1, v)
    Ik = i_k1(v)
    Iy = If(y, gfNa, ENa, EK, v)

    Ii = Is + In + Ix + Ik + Iy + I

    retcurrent = [Is, In, Ix, Ik, Iy, Ii]

    return retcurrent

def pulsefun(state0, stimpoint, period, cycletime, biasamp, stimamp):
    stimtime = stimpoint + period * cycletime
    stimmax = stimtime + 50

    t1 = np.arange(0, stimtime, 0.2)
    t2 = np.arange(stimtime, stimmax, 0.2)
    t3 = np.arange(stimmax, 5000, 0.2)

    ###Parameters

    parameters = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, biasamp]
    param = parameters[:]
    param.append(0)
    p2 = parameters[:]
    p2.append(stimamp)

    state1 = odeint(PacemakerODE, state0, t1, args=(param,), hmax=0.2)
    state2 = odeint(PacemakerODE, state1[-1, :], t2, args=(p2,), hmax=0.2)
    state3 = odeint(PacemakerODE, state2[-1, :], t3, args=(param,), hmax=0.2)

    truestate = np.concatenate([state1, state2, state3])
    truetime = np.concatenate([t1, t2, t3])

    return truestate, truetime

#parameters
parametersdep = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, -2, 0]
parametershyp = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, 0.7, 0]
parameterstimdep = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, -2, -2]
state0 = [0, 0, 0, 0, 0, 0, 0]

#Time
tmax = 10000
dt = 0.2
t = np.arange(0, tmax, dt)

#values
##Depressed model values
state = odeint(PacemakerODE, state0, t, args=(parametersdep,))
pointdiff = argrelextrema(np.array(state[:,0]), np.greater)[0] * 0.2
per = np.diff(pointdiff)[2]
period = pointdiff[3]

##Hyperpolarized values
statehyp = odeint(PacemakerODE, state0, t, args=(parametershyp,))
pointdiffhyp = argrelextrema(np.array(statehyp[:,0]), np.greater)[0] * 0.2
perhyp = np.diff(pointdiffhyp)[2]
periodhyp = pointdiffhyp[3]