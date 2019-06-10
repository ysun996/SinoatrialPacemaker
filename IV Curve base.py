from commonIV import *

currentarray = []

for v in vclamp:
    state0 = [v, 0, 0, 0, 0, 0, 0]
    state = odeint(PacemakerClamp, state0, t, args=(parametersbase,))
    currentarray.append(memcurrent(state[-1], parametersbase))

plt.plot(vclamp, currentarray)
plt.xlabel('Voltage (mV)')
plt.ylabel('Current (uA/cm^2)')
plt.title('Control')

