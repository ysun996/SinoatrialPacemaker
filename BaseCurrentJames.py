from commonJames import *

############################################################
# Currents
CurrentYN = 1
if CurrentYN == 1:
    parameters[9] = 0
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    dstate = [[0, 0, 0, 0, 0, 0]]
    for w in range(len(state[:, 0])):
        d1state = [currentVals(state[w, :], t, parameters)]
        dstate = np.concatenate((dstate, d1state), axis=0)
    Ilabels = ['Isi', 'Ina', 'Ix1', 'Ik1', 'If', 'Ii']
    for w in range(0, 5):
        plt.plot(t, dstate[1:, w], label=Ilabels[w])

    plt.legend()
    plt.ylabel('Current')
    plt.xlabel('t')
    plt.xlim(250, 1250)

############################################################
