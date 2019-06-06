
from commonJames import *
############################################################
# Action potential
ActionPotentialYN = 0
if ActionPotentialYN == 1:
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    plt.plot(t, state[:, 0])
    plt.ylabel('V')
    plt.xlabel('t')
    plt.xlim(250, 1250)
#############################################################