from commonIV import *

currentarray = []

condition = 2
for v in vclamp:
    state0 = [v, 0, 0, 0, 0, 0, 0]
    state = odeint(PacemakerODE, state0, t, args=(parametersbase, condition))
    currentarray.append(memcurrent(state[-1], parametersbase))

plt.plot(vclamp, currentarray)
plt.xlabel('Voltage (mV)')
plt.ylabel('Current (uA/cm^2)')
plt.title('Control')