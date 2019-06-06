import numpy as np
import matplotlib.pyplot as plt
import bisect
from scipy.integrate import odeint
#####################################################################################################
# This Input file was set to create a graph of individual currents.
#####################################################################################################

#####Gating variable equations and individual current equations
def dzdt(a,b,z):
    hill=a/(a+b)
    timeC=1/(a+b)
    rVal=(hill-z)/timeC
    return rVal
def am(Vm):
    max_a_m = 1
    a_m = max_a_m * ((Vm + 47) / (1 - np.exp(-(Vm + 47) / 10)))

    return a_m
def an(Vm):
    a_n=0.01*(Vm+55)/(1-np.exp(-(Vm+55/10)))
    return a_n

def As(v):
    a_s=0.05*(v+28)/(1-np.exp(-(v+28)/5))
    return a_s

def bs(v):
    b_s=0.00025*np.exp(-(v+28)/14.93)
    return b_s

def bn(Vm):
    bn=0.125*np.exp((-Vm+65)/80)
    return bn

def bm(Vm):
    max_b_m = 40
    b_m = max_b_m * np.exp(-(Vm + 72) / 17.86)

    return b_m

def ah(Vm):
    max_a_h = 0.0085
    a_h = max_a_h * np.exp(-(Vm + 71) / 5.43)

    return a_h

def bh(Vm):
    max_b_h = 2.5
    b_h = max_b_h / (1 + np.exp(-(Vm + 10) / 12.2))

    return b_h


def ad(Vm):
    max_a_d = 0.005
    a_d = max_a_d * ((Vm + 34) / (1 - np.exp(-(Vm + 34) / 10)))

    return a_d

def bd(Vm):
    max_b_d = 0.05
    b_d = max_b_d * (np.exp(-(Vm + 34) / 6.67))

    return b_d

def not_d(Vm):
    ### d' variable needed for the slow inward current
    not_d = 1/(1+np.exp(-(Vm+15)/6.67))

    return not_d

def af(Vm):
    max_a_f = 0.002468
    a_f = max_a_f * np.exp(-(Vm + 47) / 20)

    return a_f

def bf(Vm):
    max_b_f = 0.05
    b_f = max_b_f / (1 + np.exp(-(Vm + 13) / 11.49))

    return b_f


def ax1(Vm):
    max_a_x1 = 0.0025
    a_x1 = (max_a_x1 * np.exp((Vm + 30) / 12.11)) / (1 + np.exp((Vm + 30) / 50))

    return a_x1

def bx1(Vm):
    max_b_x1 = 0.0065
    b_x1 = max_b_x1 * (np.exp(-(Vm - 20) / 16.67)) / (1 + np.exp(-(Vm - 20) / 25))

    return b_x1

def l_bar_x1(Vm):
    l_bar = 2.25 * ((np.exp((Vm + 95)/25)-1)/(np.exp((Vm + 45)/25)))

    return l_bar




def ay(Vm):
    #Polynomial Regression
    w=0.03375000000000002



    if isinstance(Vm,float)==True:
        a = [1.4557319134e-12,4.0945641782e-10,4.6549818992e-08,2.4903140216e-06,6.1460577425e-05,4.7453248494e-04,\
             2.5019715465e-03]
        a_y = a[0]*Vm**6+a[1]*Vm**5+a[2]*Vm**4+a[3]*Vm**3+a[4]*Vm**2+a[5]*Vm+a[6]
        if a_y <= 0:
            a_y = 0.00001
    else:
        a = [1.4557319134e-12,4.0945641782e-10,4.6549818992e-08,2.4903140216e-06,6.1460577425e-05,4.7453248494e-04,\
             2.5019715465e-03]
        a_y = a[0]*Vm**6+a[1]*Vm**5+a[2]*Vm**4+a[3]*Vm**3+a[4]*Vm**2+a[5]*Vm+a[6]
        for i in range(len(Vm)):
            if a_y[i] <= 0:
                a_y[i] = 0.00001        
    

    a_y = a_y*w
    #Modified version of the original equation
    #max_a_y = 0.00005
    #a_y = (max_a_y * np.exp(-(Vm + 14.25) / 14.93))/3
    return a_y

def by(Vm):

    #Polynomial Regression
    w=0.03375000000000002
    
    if isinstance(Vm,float)==True:

        a = [3.5607174324e-13,3.9587887660e-11,-6.9345321240e-09,-8.8541673551e-07,4.5605591007e-05,9.4190808268e-03,\
             3.3771510156e-01]
        b_y = a[0]*Vm**6+a[1]*Vm**5+a[2]*Vm**4+a[3]*Vm**3+a[4]*Vm**2+a[5]*Vm+a[6]
        if b_y <= 0:
            b_y = 0.00001
    else:
        a = [3.5607174324e-13,3.9587887660e-11,-6.9345321240e-09,-8.8541673551e-07,4.5605591007e-05,9.4190808268e-03,\
            3.3771510156e-01]
        b_y = a[0]*Vm**6+a[1]*Vm**5+a[2]*Vm**4+a[3]*Vm**3+a[4]*Vm**2+a[5]*Vm+a[6]
        for i in range(len(Vm)):
            if b_y[i] <= 0:
                b_y[i] = 0.00001
    b_y = b_y*w    
    
    #Modified version of the original equation
    #max_b_y = 0.001
    #b_y = (max_b_y * (Vm + 14.25) / (1 - np.exp(-(Vm + 14.25) / 5)))/3
    
    return b_y

#####Time-independent background potassium current (I_k1)

def i_k1(Vm):
    IPot=1.3*(np.exp((Vm+110)/25)-1)/(np.exp((Vm+60)/12.5)+np.exp((Vm+60)/25))+\
         0.26*(Vm+30)/(1-np.exp(-(Vm+30)/25))
    return IPot

#Parameters
gNa = 4.4
gnNa = 0.066
gSi = 0.5175
gnSi = 0.161
gfNa = 1.2
ENa = 40
ESi = 70
EK = -93
Cm = 6
I = 0
parameters = [gNa, gnNa, gSi, gnSi, gfNa, ENa, ESi, EK, Cm, I]



###Pacemaker Function

def PacemakerODE(state, t, parameters):
    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    ENa = parameters[5]
    ESi = parameters[6]
    EK = parameters[7]
    Cm = parameters[8]
    I = parameters[9]



    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]



    Isi = (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m**3)*h + gnNa) * (v - ENa)
    Ix1 = x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = (y**2 * gfNa) * (v - ENa) + (v-EK) * y**2 * (-120/EK)
    Ii = Isi + Ina + Ix1 + Ik1 + If + I 
    dv = -Ii/Cm


    dd = (ad(v) * (1-d)) - (bd(v) * d)
    df = (af(v) * (1-f)) - (bf(v) * f)
    dm = (am(v) * (1-m)) - (bm(v) * m)
    dh = (ah(v) * (1-h)) - (bh(v) * h)
    dx1 = (ax1(v) * (1-x1)) - (bx1(v) * x1)
    dy = (ay(v) * (1-y)) - (by(v) * y)
    rstate = [dv, dd, df, dm, dh, dx1, dy]
    return rstate

#Function to find individual current values of PacemakerODE
def currentVals(state, t, parameters):
    gNa = parameters[0]
    gnNa = parameters[1]
    gSi = parameters[2]
    gnSi = parameters[3]
    gfNa = parameters[4]
    ENa = parameters[5]
    ESi = parameters[6]
    EK = parameters[7]
    Cm = parameters[8]
    I = parameters[9]


    v = state[0]
    d = state[1]
    f = state[2]
    m = state[3]
    h = state[4]
    x1 = state[5]
    y = state[6]


    Isi = (gSi * d * f + gnSi * not_d(v)) * (v - ESi)
    Ina = (gNa * (m**3)*h + gnNa) * (v - ENa)
    Ix1 = x1 * l_bar_x1(v)
    Ik1 = i_k1(v)
    If = y**2*(gfNa * (v - ENa) + (v-EK) * (-120/EK))
    Ii = Isi + Ina + Ix1 + Ik1 + If + I 
    return [Isi,Ina,Ix1,Ik1,If,Ii]

#Function to find local maxima and cycle length
def cycle_len(t,s):
    logic1 = True
    localMax = [0]
    for i in range(2500,5000):
        if logic1 and s[i,0] > -20 and s[i,0] < s[(i-1),0]:
            localMax.append(i - 1)
            logic1=False
        elif s[i,0] > -20 and s[(i-1),0] < -20:
            logic1=True
    cyclelen = t[localMax[2]]-t[localMax[1]]
    return cyclelen,localMax

#Pulse function
def pulseFnc(s0,t,p,phase,stim):
    for x in [0,1]:
        if x == 0:
            tIndexL=bisect.bisect_left(t,phase)
            t1=t[0:tIndexL]
            tplot=t1
            s = odeint(PacemakerODE, s0, t1, args=(parameters,))
            stateFinal = s
        else:
            CycleValue=len(t)
            stimCycleValue=phase + 50
            tIndexR=bisect.bisect_left(t,CycleValue)
            tIndexStim=bisect.bisect_left(t,stimCycleValue)
            p[9] = stim
            t1=t[tIndexL:tIndexStim]
            tplot=np.append(tplot,t1)
            state1=s[-1,:]
            s = odeint(PacemakerODE, state1, t1, args=(parameters,))
            stateFinal = np.concatenate((stateFinal,s),axis=0)
            p[9]=0
            t1=t[tIndexStim:tIndexR]
            tplot=np.append(tplot,t1)
            state1=s[-1,:]
            s = odeint(PacemakerODE, state1, t1, args=(parameters,))
            stateFinal = np.concatenate((stateFinal,s),axis=0)
            tIndexL=tIndexR
    return stateFinal
            
#Function to find Time Constant
def timeC(a,b):
    r = 1/(a+b)
    return r

#Initial state
state0=[-1.26312652e+01,  6.73685294e-01,  5.28447665e-01,  9.60815694e-01,
  4.83891944e-07,  3.70080101e-02,  1.24208141e-01]
           
parameters[9] = 0

#Time
tmax = 2019
dt = 0.4
t = np.arange(0, tmax, dt)


#########################################################################################
# Phase Response Curves. Set HyperPolarYN to 1 for hyperpolarized phase-response curve or 0 for depolarized phase response
# Curve.
PhaseResetYN = 0
HyperPolarYN = 1
if HyperPolarYN == 1:
    PolPulse = [0.5,2.0,3.5]
elif HyperPolarYN == 0:
    PolPulse = [-0.6,-1.3,-1.9,-2.8]
if PhaseResetYN == 1:
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    cycleLength, localMx = cycle_len(t,state)
    print(cycleLength)
    phasePercent = np.arange(0.05,1,0.025)
    for w in PolPulse:
        print('now on '+str(w))
        phaseChange = np.zeros(len(phasePercent))
        for i in range(len(phasePercent)):
            phase = t[localMx[1]]+cycleLength*phasePercent[i]
            statePulsed = pulseFnc(state0,t,parameters,phase,w)
            cyclePulsedLength, localPulsedMx = cycle_len(t,statePulsed)
            phaseChange[i] = (cyclePulsedLength-cycleLength)/cycleLength*100
            
        plt.plot(phasePercent*100,phaseChange, 'o',label=str(w) + '\u03BC'+'A/c$m^2$')
    plt.legend()
    plt.ylabel('\u0394'+'$\phi$'+' (%)')
    plt.xlabel('% of Cycle ($\phi$)')
#############################################################################################

#############################################################################################
#v vs. dv
vdvYN = 0
if vdvYN == 1:
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    dstate=[[0,0,0,0,0,0,0]]
    for w in range(len(state[:,0])):
        d1state=[PacemakerODE(state[w,:],t,parameters)]
        dstate=np.concatenate((dstate,d1state),axis=0)
    plt.plot(state[:,0],dstate[1:,0])
    plt.ylabel('dV/dt')
    plt.xlabel('V')

#############################################################################################   

############################################################
# Currents
CurrentYN = 1
if CurrentYN == 1:
    parameters[9] = 0
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    dstate=[[0,0,0,0,0,0]]
    for w in range(len(state[:,0])):
        d1state=[currentVals(state[w,:],t,parameters)]
        dstate=np.concatenate((dstate,d1state),axis=0)
    Ilabels=['Isi','Ina','Ix1','Ik1','If','Ii']
    for w in range(0,5):
        plt.plot(t,dstate[1:,w],label=Ilabels[w])
    
    plt.legend()
    plt.ylabel('Current')
    plt.xlabel('t')
    plt.xlim(250,1250)


############################################################


############################################################
#Action potential
ActionPotentialYN = 0
if ActionPotentialYN == 1:
    state = odeint(PacemakerODE, state0, t, args=(parameters,))
    plt.plot(t,state[:,0])
    plt.ylabel('V')
    plt.xlabel('t')
    plt.xlim(250,1250)
#############################################################