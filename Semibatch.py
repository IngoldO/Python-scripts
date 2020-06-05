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
def model(y,t,F,caF,rho,cp,dHr,Tf,k0,Ea,Tstart,Tset,Kc,Ki,Qmax):
        Na = y[0]    # Concentration of component A (mol/L)
        Nb = y[1]    # Concentration of component B (mol/L)
      
        V = y[3]     # Reactor volume (L)
        Tfun = y[4]     # Reactor temperature (K)
        eint = y[5]     # Integrated control error (K*s)
        T_error = Tfun-Tset

        #if Tfun <= Tset:
        #        cont = abs(Kc*T_error+Ki*eint) # PI controller value
        #else:
        cont = Kc*T_error+Ki*eint 

        if cont==0:
                Q = 0
        else: 
                Q = Qmax*cont/abs(1+cont) # Heating / cooling power (W)
        #else:
         #       Q = Qmax*cont/(-1+cont)
        #print("Temp = " + str(Tfun))
        #print("cont = " + str(cont))
        #print("Q = " + str(Q))
        

        if Tfun < Tstart:
                F = 0   # Feed only starts at starting temperature
        k = k0*math.exp(-Ea/(8.3145*Tfun))        # Reaction rate constant (L/mol/s)
        dNadt =  F*caF-k*Na*Nb/V     # molar amount change of component A [mol/L/s]
        dNbdt =  -k*Na*Nb/V          # molar amount change of component A [mol/L/s]
        dNcdt =  k*Na*Nb/V           # molar amount change of component A [mol/L/s]
        dVdt =  F                    # Volume change rate (L/s)
        dTdt = (F*rho*cp*(Tf-Tfun)+k*Na*Nb*dHr/V+Q)/(V*rho*cp) # rate of temperature change [K]
        deintdt = T_error              # Integration of control error (K*s)
        #print(k)
        return [dNadt,dNbdt,dNcdt,dVdt,dTdt,deintdt]
 
# Parameters
F = 0.5     # Feed rate (L/s)
caF = 5 # Feed concentration of A (L/s)
Tf = 298  # Temperature of the feed (K)
Tstart = 298  # Starting temperature of the feed [K]
Tset = 273+80 # Setpoint temperature of PI controller
Kc = -1 # PI controller proportional constant
Ki = -0.0001 # PI controller integral constant
Qmax = 1000000 # Maximal heating / cooling power (W)

# intial conditions
ca0 = 0.0    # Concentration of component A (mol/L)
cb0 = 2    # Concentration of component B (mol/L)
cc0 = 0      # Concentration of component C (mol/L)
V0 = 2000    # Reactor volume (L)
T0 = 298     # Starting temperature (K)
k0 = 10**9   # Arrhenius constant (L/mol/s)
Ea = 13000    # activation energy (J/mol)
p0 = Kc*(T0-Tset)   # Initial heating/cooling (W)
y0 = [ca0*V0,cb0*V0,cc0*V0,V0,T0,p0] # vector of initial conditions

#Physical data
rho = 1       # density of the solvent (kg/L)
dHr = 140000 # enthalpy of reaction (J/mol)
cp = 4186     # heat capacity of solvent (J/kg/K)



# time points
t = np.linspace(0,1800)

# solving ODE
y = odeint(model,y0,t,args=(F,caF,rho,cp,dHr,Tf,k0,Ea,Tstart,Tset,Kc,Ki,Qmax))

Na = y[:,0]
Nb = y[:,1]
Nc = y[:,2]
V = y[:,3]
T = y[:,4]
p = y[:,5]

T_error = T-Tset
cont = Kc*T_error+Ki*np.trapz(T_error,t)
Q = Qmax*cont/abs(1+cont)

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
plt.plot(t,T-273.15)
plt.xlabel('Time (s)')
plt.ylabel('Temperature (deg C)')

plt.figure()
plt.plot(t,T-273.15)
plt.plot(t,p)
plt.xlabel('Time (s)')
plt.ylabel('Integrated temperature error (deg C)')
