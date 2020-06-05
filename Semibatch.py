# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 18:32:13 2017
 
@author: ingoldo
"""
 
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
 
# function returning dydt
def model(y,t,F,caF,rho,cp,dHr,Q,Tf,k0,Ea,Tstart,Tset,Kc,Ki):
        Na = y[0]    # Concentration of component A (mol/L)
        Nb = y[1]    # Concentration of component B (mol/L)
        V = y[3]     # Reactor volume (L)
        T = y[4]     # Reactor temperature (K)
        Q = y[5]     # Heating / cooling power (W)
        if T < Tstart:
                F = 0   # Feed only starts at starting temperature
        k = k0*math.exp(-Ea/(8.3145*T))        # Reaction rate constant (L/mol/s)
        dNadt =  F*caF-k*Na*Nb/V     # molar amount change of component A [mol/L/s]
        dNbdt =  -k*Na*Nb/V          # molar amount change of component A [mol/L/s]
        dNcdt =  k*Na*Nb/V           # molar amount change of component A [mol/L/s]
        dVdt =  F                    # Volume change rate (L/s)
        dTdt = (F*rho*cp*(Tf-T)+k*Na*Nb*dHr/V+Q)/(V*rho*cp) # rate of temperature change [K]
        dQdt = Kc+Ki*(T-Tset)
 #       dQdt = 0.0001
        return [dNadt,dNbdt,dNcdt,dVdt,dTdt,dQdt]
 
# intial conditions
ca0 = 0.0    # Concentration of component A (mol/L)
cb0 = 0.1    # Concentration of component B (mol/L)
cc0 = 0      # Concentration of component C (mol/L)
V0 = 6000    # Reactor volume (L)
T0 = 298     # Starting temperature (K)
k0 = 1        # Arrhenius constant (L/mol/s)
Ea = 1.6*10**5 # activation energy (J/mol)
Q0 = 0       # Initial heating/cooling (W)
y0 = [ca0*V0,cb0*V0,cc0*V0,V0,T0,Q0] # vector of initial conditions

#Physical data
rho = 1       # density of the solvent (kg/L)
dHr = 1400000 # enthalpy of reaction (J/mol)
cp = 4186     # heat capacity of solvent (J/kg/K)

# Parameters
F = 5     # Feed rate (L/s)
caF = 0.1 # Feed concentration of A (L/s)
Tf = 298  # Temperature of the feed (K)
Q = 20000 # cooling power (W)
Tstart = 298  # Starting temperature of the feed [K]
Tset = 273+80 # Setpoint temperature of PI controller
Kc = -1 # PI controller proportional constant
Ki = -1 # PI controller integral constant


# time points
t = np.linspace(0,3900)

# solving ODE
y = odeint(model,y0,t,args=(F,caF,rho,cp,dHr,Q,Tf,k0,Ea,Tstart,Tset,Kc,Ki))

Na = y[:,0]
Nb = y[:,1]
Nc = y[:,2]
V = y[:,3]
T = y[:,4]
ca = Na/V
cb = Nb/V
cc = Nc/V

 
# plot results
plt.close
plt.plot(t,ca)
plt.plot(t,cb)
plt.plot(t,cc)
plt.xlabel('Time (s)')
plt.ylabel('Concentration (g/L)')
plt.legend(['c_A','c_B','c_C'])
plt.show
plt.close

plt.figure()
plt.plot(t,T)