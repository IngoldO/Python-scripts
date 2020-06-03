# -*- coding: utf-8 -*-
"""
Created on Sat Dec 16 18:32:13 2017
 
@author: ingoldo
"""
 
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
 
# function returning dydt
def model(y,t,F,Sf):
        mumax = 0.3 # maximal growth rate (1/s)
        YXs = 0.1     # yield coefficient
        X = y[0]    # cell concentration (g/L)
        S = y[1]    # Substrate concentration (g/L)
        Ks = 1      # Half saturation constant (g/L)
        b = 0.2     # decay rate coefficient (1/s)
        V = 100     # Bioreactor volume (L)
        mu = mumax*S/(Ks+S) # reaction rate constant (substrate inhibition)
        dXdt = mu*X-b*X     # growth rate of cells (g/L/h)
        dSdt = -YXs*X + F * Sf / V # Substrate change rate (g/L/h)
        return [dXdt,dSdt]
 
# intial condition
y0 = 5
 
# time points
t = np.linspace(0,80)
 
# solving ODE
y = odeint(model,[0.1,10],t,args=(0.1,0.1,))

 
# plot results
plt.close
plt.plot(t,y)
plt.xlabel('Time (h)')
plt.ylabel('Concentration (g/L)')
plt.show
plt.close