from commonDep import *

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

state0 = [0, 0, 0, 0, 0, 0, 0]
state = odeint(PacemakerODE, state0, t, args=(parameters,))
pointdiff = argrelextrema(np.array(state[:,0]), np.greater)[0] * 0.2
per = np.diff(pointdiff)

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

stim = pulsefun(period, per, 0.29, -2, -2)

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