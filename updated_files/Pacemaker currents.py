from updated_files.common import *

condition = 0

state0 = [-66, 0, 0, 0, 0, 0, 0]

state = odeint(PacemakerODE, state0, t, args=(parametersbase,condition))
dv = []
for points in state:
    pacemaker = PacemakerODE(points, t, parametersbase, condition)
    dv.append(pacemaker[0])
plt.plot(t, state[:,0])
plt.figure()
plt.plot(state[:,0], dv)









