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

    Isi = 0.78 * Isi(gSi, gnSi, d, f, v, ESi)
    Ina = Ina(gNa, m, h, gnNa, v, ENa)
    Ix1 = 0.8 * Ix1(x1, v)
    Ik1 = i_k1(v)
    If = If(y, gfNa, ENa, EK, v)

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

    Isi = 0.78 * Isi(gSi, gnSi, d, f, v, ESi)
    Ina = Ina(gNa, m, h, gnNa, v, ENa)
    Ix1 = 0.8 * Ix1(x1, v)
    Ik1 = i_k1(v)
    If = If(y, gfNa, ENa, EK, v)

    Ii = Isi + Ina + Ix1 + Ik1 + If + I

    retcurrent = [Isi, Ina, Ix1, Ik1, If, Ii]

    return retcurrent

def pulsefun(stimpoint, period, cycletime, biasamp, stimamp):
    stimtime = stimpoint + period * cycletime
    stimmax = stimtime + 50

    t1 = np.arange(0, stimtime, 0.2)
    t2 = np.arange(stimtime, stimmax, 0.2)
    t3 = np.arange(stimmax, 5000, 0.2)

    ###Parameters

    parameters = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, biasamp]
    param = parameters[:]
    p1.append(0)
    p2 = parameters[:]
    p2.append(stimamp)


    state1 = odeint(PacemakerODE, state0, t1, args=(param,), hmax=0.2)
    state2 = odeint(PacemakerODE, state1[-1, :], t2, args=(p2,), hmax=0.2)
    state3 = odeint(PacemakerODE, state2[-1, :], t3, args=(param,), hmax=0.2)

    truestate = np.concatenate([state1[:, 0], state2[:, 0], state3[:, 0]])
    truetime = np.concatenate([t1, t2, t3])

    plotfun = [truestate, truetime]

    return plotfun

#parameters
parametersdep = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, -2, 0]
state0 = [0, 0, 0, 0, 0, 0, 0]

#Time
tmax = 10000
dt = 0.2
t = np.arange(0, tmax, dt)

#values
state = odeint(PacemakerODE, state0, t, args=(parametersdep,))
pointdiff = argrelextrema(np.array(state[:,0]), np.greater)[0] * 0.2
per = np.diff(pointdiff)
period = pointdiff[3]