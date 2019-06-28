from updated_files.commonJames import *


# Initial state


#########################################################################################
# Phase Response Curves. Set HyperPolarYN to 1 for hyperpolarized phase-response curve or 0 for depolarized phase response
# Curve.
PhaseResetYN = 1
HyperPolarYN = 1
if HyperPolarYN == 1:
    PolPulse = [0.5, 2.0, 3.5]
elif HyperPolarYN == 0:
    PolPulse = [-0.6, -1.3, -1.9, -2.8]
if PhaseResetYN == 1:
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    cycleLength, localMx = cycle_len(t, state)
    print(cycleLength)
    phasePercent = np.arange(0.05, 1, 0.025)
    for w in PolPulse:
        print('now on ' + str(w))
        phaseChange = np.zeros(len(phasePercent))
        for i in range(len(phasePercent)):
            phase = t[localMx[1]] + cycleLength * phasePercent[i]
            statePulsed = pulseFnc(state0, t, parameters, phase, w)
            cyclePulsedLength, localPulsedMx = cycle_len(t, statePulsed)
            phaseChange[i] = (cyclePulsedLength - cycleLength) / cycleLength * 100

        plt.plot(phasePercent * 100, phaseChange, 'o', label=str(w) + '\u03BC' + 'A/c$m^2$')
    plt.legend()
    plt.ylabel('\u0394' + '$\phi$' + ' (%)')
    plt.xlabel('% of Cycle ($\phi$)')
#############################################################################################

