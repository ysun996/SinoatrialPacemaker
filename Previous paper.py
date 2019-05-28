import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
aMax=[1,0.0085,0.005,0.002468,0.05,0.00025]
bMax=[40,2.5,0.05,0.05,0.0025,0.0065]
#m, h, d, f, s, x1
Gmax=[5,0.075,0.45,0.14]
#Na,Na',Si,Si'
RevPotential=[40,70]
#Na,Si
def am(vm,aMax):
    a_m=(aMax*(vm+47))/(1-np.exp(-(vm+47)/10))
    return a_m

def bm(vm,bMax):
    b_m=bMax*np.exp(-(vm+72)/17.86)
    return b_m

def ah(vm,aMax):
    a_h=aMax*np.exp(-(vm+71)/5.43)
    return a_h

def bh(vm,bMax):
    b_h=bMax/(1+np.exp(-(vm+10)/12))
    return b_h

def ad(vm,aMax):
    a_d=(aMax*(vm+34))/(1-np.exp(-(vm+34)/10))
    return a_d

def bd(vm,bMax):
    b_d=bMax*np.exp(-(vm+34)/6.67)
    return b_d

def af(vm,aMax):
    a_f=aMax*np.exp(-(vm+47)/20)
    return a_f

def bf(vm,bMax):
    b_f=bMax/(1+np.exp(-(vm+13)/11.49))
    return b_f

def As(vm,aMax):
    a_s=(aMax*(vm+28))/(1-np.exp(-(vm+28)/5))
    return a_s

def bs(vm,bMax):
    b_s=bMax*np.exp(-(vm+28)/14.93)
    return b_s

def ax1(vm,Amax):
    a_x1=Amax*np.exp((vm+30/12.11))/(1+np.exp((vm+30)/50))
    return a_x1

def bx1(vm,Bmax):
    b_x1=Bmax*np.exp(-(vm-20/16.67))/(1+np.exp(-(vm-20)/25))
    return b_x1

def ik2(vm):
    i_k2=0.6*(np.exp((vm+110)/25) -1)/(np.exp((vm+60)/12.5)+np.exp((vm+60)/25))
    return i_k2

def ix1(vm):
    i_x1=2.25*(np.exp((vm+95)/25) -1)/np.exp((vm+45)/25)
    return i_x1

def notD(vm):
    dnot=1/(1+np.exp(-(vm+15)/6.67))
    return dnot

def dzdt(z,a,b):
    dz=a*(1-z)+b*z
    return dz

def gNa(m,h,Gmax):
    gNaOut=Gmax[0]*(m**3)*h+Gmax[1]
    return gNaOut

def gSi(d,f,nd,Gmax):
    gSiOut=Gmax[2]*d*f+Gmax[3]*d
    return gSiOut

def Ina(G_Na,vm,E):
    InOut=G_Na*(vm-E[0])
    return InOut

def Isi(G_Si,vm,E):
    IsOut=G_Si*(vm-E[1])
    return IsOut

def IK2(s,ik):
    I_K2=s*ik
    return I_K2

def IX1(x,ix):
    I_X1=x*ix
    return I_X1

def IK1(vm,ik):
    I_K1=ik/0.6+((vm+30)/5)/(1-np.exp((vm+30)/25))
    return I_K1

def PacemakerODE(state,t,parameters):
    aMax=parameters[0]
    bMax=parameters[1]
    gMax=parameters[2]
    e=parameters[3]
    Cm=6
    v=state[0]
    m=state[1]
    h=state[2]
    d=state[3]
    f=state[4]
    s=state[5]
    x1=state[6]
    I_i=Ina(gNa(m,h,gMax),v,e)+Isi(gSi(d,f,notD(v),gMax),v,e)+IK2(s,ik2(v))+IX1(x1,ix1(v))+IK1(v,ik2(v))
    dv=-I_i/Cm
    dm=dzdt(m,am(v,aMax[0]),bm(v,bMax[0]))
    dh=dzdt(h,ah(v,aMax[1]),bh(v,bMax[1]))
    dd=dzdt(d,ad(v,aMax[2]),bd(v,bMax[2]))
    df=dzdt(f,af(v,aMax[3]),bf(v,bMax[3]))
    ds=dzdt(s,As(v,aMax[4]),bs(v,bMax[4]))
    dx1=dzdt(x1,ax1(v,aMax[5]),bx1(v,bMax[5]))
    rstate=[dv,dm,dh,dd,df,ds,dx1]
    return rstate
#Time
tmax = 1000
dt = 0.2
t = np.arange(0, tmax, dt)

# m = np.arange(0, 1, 0.005)
# # h = 0
# v = np.arange(-100, 100, 1)
# # dm=dzdt(m,am(voltage, Amax[0]),bm(voltage,Bmax[0]))
# # dh=dzdt(h,ah(voltage,Amax[1]),bh(voltage,Bmax[1]))
#
# derivarray = [dzdt(m,am(v,aMax[0]),bm(v,bMax[0])),dzdt(m,ah(v,aMax[1]),bh(v,bMax[1])), dzdt(m,ad(v,aMax[2]),bd(v,bMax[2])),
#              dzdt(m,af(v,aMax[3]),bf(v,bMax[3])), dzdt(m,As(v,aMax[4]),bs(v,bMax[4])), dzdt(m,ax1(v,aMax[5]),bx1(v,bMax[5]))]
#
# for i in range(len(derivarray)):
#     # plt.subplot(len(derivarray), 1, i+1)
#     for j in v:
#         volts = np.full(len(derivarray[i]), j)
#         plt.plot(volts, derivarray[i])
#     plt.figure()



parameters=[aMax,bMax,Gmax,RevPotential]
state0=[40,0,0.3,0,0,0,0]

state=odeint(PacemakerODE,state0,t,args=(parameters,))
plt.plot(t,state[:,0])