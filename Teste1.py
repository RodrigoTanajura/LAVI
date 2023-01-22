#%%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

w = np.linspace(0,50,1001)
k = 10000
c = 100
u = 10
m = 10
m0 = 1
y0 = [0,0]
t = 1
wn = ((k/m)**(1/2))
r = w/wn

def integrax(y, w, k, c, m, u):
    x, x_dot = y
    dydt = [x_dot,(m*u*w**2 *np.cos(w*t)-(k*x)-(c*x_dot))/m]
    return dydt

def integray(y, w, k, c, m, u):
    x, x_dot = y
    dydt = [x_dot,(m*u*w**2 *np.sin(w*t)-(k*x)-(c*x_dot))/m]
    return dydt

solx = odeint(integrax, y0, w, args=(k, c, m, u))
soly = odeint(integray, y0, w, args=(k, c, m, u))

plt.plot(r,solx[:,0])
plt.xlabel('w/wn')
plt.ylabel('Amplitude')
plt.show()