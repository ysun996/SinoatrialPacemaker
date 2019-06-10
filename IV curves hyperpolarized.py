from commonIV import *

currentarrayh = []

for v in vclamp:
    state0 = [v, 0, 0, 0, 0, 0, 0]
    state = odeint(PacemakerClamp, state0, t, args=(parametershyp,))
    currentarrayh.append(memcurrent(state[-1], parametershyp))

plt.plot(vclamp, currentarrayh, '-')
plt.xlabel('Voltage (mV)')
plt.ylabel('Current (uA/cm^2)')
plt.title('Hyperpolarized model')