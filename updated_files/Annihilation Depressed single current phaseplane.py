from updated_files.commonDep import *


dv = []
stimpulse, timepulse = pulsefun(state0, period, per, 0.29, -2, -2)
stimtime = period + per * 0.29

for points in range(len(timepulse)):
    if timepulse[points] < stimtime:
        pacemaker = PacemakerODE(stimpulse[points], t, parametersdep, condition)
    elif timepulse[points] >= stimtime and timepulse[points] <= stimtime + 50:
        pacemaker = PacemakerODE(stimpulse[points], t, parameterstimdep, condition)
    else:
        pacemaker = PacemakerODE(stimpulse[points], t, parametersdep, condition)
    dv.append(pacemaker[0])

plt.plot(stimpulse[:, 0], dv)
plt.title('Phase-plane relationship at 29%')