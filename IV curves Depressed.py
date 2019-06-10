from commonIV import *

currentarrayd = []

for v in vclamp:
    state0 = [v, 0, 0, 0, 0, 0, 0]
    state = odeint(PacemakerClamp, state0, t, args=(parametersdep,))
    currentarrayd.append(memcurrent(state[-1], parametersdep))

plt.plot(vclamp, currentarrayd)
plt.xlabel('Voltage (mV)')
plt.ylabel('Current (uA/cm^2)')
plt.title('Depressed model')