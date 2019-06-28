from updated_files.common import *
from scipy.integrate import odeint
from scipy.signal import argrelextrema


def pulsefun(state0, stimpoint, period, cycletime, biasamp, stimamp, condition = 0):
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

    state1 = odeint(PacemakerODE, state0, t1, args=(param, condition), hmax=0.2)
    state2 = odeint(PacemakerODE, state1[-1, :], t2, args=(p2, condition), hmax=0.2)
    state3 = odeint(PacemakerODE, state2[-1, :], t3, args=(param, condition), hmax=0.2)

    truestate = np.concatenate([state1, state2, state3])
    truetime = np.concatenate([t1, t2, t3])

    return truestate, truetime

#parameters
parameterstimdep = [4.4, 0.066, 0.5175, 0.161, 1.2, 40, 70, -93, 6, -2, -2]
state0 = [0, 0, 0, 0, 0, 0, 0]


#values
##PacemakerODE condition
condition = 0

##Depressed model values
state = odeint(PacemakerODE, state0, t, args=(parametersdep, condition))
pointdiff = argrelextrema(np.array(state[:,0]), np.greater)[0] * 0.2
per = np.diff(pointdiff)[2]
period = pointdiff[3]

##Hyperpolarized values
statehyp = odeint(PacemakerODE, state0, t, args=(parametershyp, condition))
pointdiffhyp = argrelextrema(np.array(statehyp[:,0]), np.greater)[0] * 0.2
perhyp = np.diff(pointdiffhyp)[2]
periodhyp = pointdiffhyp[3]

