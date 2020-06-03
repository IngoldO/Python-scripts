# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 18:32:13 2017
 
@author: ingoldo
"""
 
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
 
# function returning dydt
def model(y,t,F,caF):
        Na = y[0]    # Concentration of component A (mol/L)
        Nb = y[1]    # Concentration of component B (mol/L)
        V = y[3]     # Reactor volume (L)
        k = 1      # Reaction rate constant (L/mol/s)
        dNadt =  F*caF-k*Na*Nb/V     # molar amount change of component A [mol/L/s]
        dNbdt =  -k*Na*Nb/V          # molar amount change of component A [mol/L/s]
        dNcdt =  k*Na*Nb/V           # molar amount change of component A [mol/L/s]
        dVdt =  F                    # Volume change rate (L/s)
        return [dNadt,dNbdt,dNcdt,dVdt]
 
# intial conditions
ca0 = 0.0    # Concentration of component A (mol/L)
cb0 = 0.1    # Concentration of component B (mol/L)
cc0 = 0      # Concentration of component C (mol/L)
V0 = 6000    # Reactor volume (L)
y0 = [ca0*V0,cb0*V0,cc0*V0,V0] # vector of initial conditions

# Parameters
F = 1   # Feed rate (L/s)
caF = 0.1 # Feed concentration of A (L/s)
 
# time points
t = np.linspace(0,1900)
 
# solving ODE
y = odeint(model,y0,t,args=(F,caF))

Na = y[:,0]
Nb = y[:,1]
Nc = y[:,2]
V = y[:,3]
ca = Na/V
cb = Nb/V
cc = Nc/V

 
# plot results
plt.close
plt.plot(t,cb)
plt.xlabel('Time (s)')
plt.ylabel('Concentration (g/L)')
plt.show
plt.close