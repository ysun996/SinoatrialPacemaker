from commonJames import *

#############################################################################################
# v vs. dv
vdvYN = 0
if vdvYN == 1:
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    dstate = [[0, 0, 0, 0, 0, 0, 0]]
    for w in range(len(state[:, 0])):
        d1state = [PacemakerODE(state[w, :], t, parameters)]
        dstate = np.concatenate((dstate, d1state), axis=0)
    plt.plot(state[:, 0], dstate[1:, 0])
    plt.ylabel('dV/dt')
    plt.xlabel('V')

#############################################################################################

