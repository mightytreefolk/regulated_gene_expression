import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import random
import math as m

n_A = 6.023E23 #Avogadro's Number

# This function models the ODE of an unregulated genetic expression of a prokaryotic cell 
# k0 is the rate at which mRNA is being generated
# k1 is the rate at which proteins are being 
# delm and delp are the degradation constants of mRNA and proteins respectively
def unregulated_gene_expression(z,t):
    k0 = 0.0001
    k1 = 0.0001
    delm = 0.03
    delp = 0.0002
    m0 = z[0]
    p0 = z[1]
    dmdt = k0 - delm*m0
    dpdt = k1*m0 - delp*p0

    dzdt = dmdt, dpdt
    return dzdt

# initial condition
z0 = [0.1, 0.3]

# number of time points
n = 10000

# time points
t = np.linspace(0,1000,n)

# store solutions
m = np.empty_like(t)
p = np.empty_like(t)

# record initial conditions
m[0] = z0[0]
p[0] = z0[1]

# solve ODE
for i in range(1,n):
    
    # span for next time step
    tspan = [t[i-1],t[i]]
    
    # solve for next step
    z = odeint(unregulated_gene_expression,z0,tspan)
    
    # store solution for plotting
    m[i] = z[1][0]
    p[i] = z[1][1]
    
    # next initial condition
    z0 = z[1] 

plt.subplot(211)
plt.plot(t,m, 'r--', label='A(t)')
plt.grid(True)
plt.legend(loc='best')

plt.subplot(212)
plt.plot(t,p, 'r--', label='R(t)')
plt.grid(True)
plt.legend(loc='best')

plt.show()

