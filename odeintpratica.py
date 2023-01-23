#%%
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

def pendx(y, t, k, c, m, u, m0, k2):
    x, x_dot = y
    dydt = [x_dot,m0*u*w**2 *np.cos(w*t)-(k*x/m)-((k2*x**3)/m)-(c*x_dot/m)]
    return dydt
def pendy(y, t, k, c, m, u, m0, k2):
    x, x_dot = y
    dydt = [x_dot,m0*u*(w**2)*np.sin(w*t)-(k*x/m)-((k2*x**3)/m)-(c*x_dot/m)]
    return dydt

t = np.linspace(0, 100, 100001)
w = 10
k = 10000
k2 = 10000
c = 100
u = 10
m = 10
m0 = 1
y0 = [1,1]

solx = odeint(pendx, y0, t, args=(k, c, m, u, m0, k2))
soly = odeint(pendy, y0, t, args=(k, c, m, u, m0, k2))
plt.plot(solx[:,0], soly[:, 0], 'b', label='x(t) vs y(t)')
plt.legend(loc='best')
plt.xlabel('x(t)')
plt.ylabel('y(t)')
plt.title('Numérica')
plt.show()
t = np.linspace(0, 3, 1001)
solx = odeint(pendx, y0, t, args=(k, c, m, u, m0, k2))
soly = odeint(pendy, y0, t, args=(k, c, m, u, m0, k2))
plt.plot(t, solx[:,0],'b', label='x(t)')
plt.plot(t, soly[:,0],'r', label='y(t)')
plt.legend(loc='best')
plt.xlabel('t')
plt.ylabel('Amplitude')
plt.title('Numérica')
plt.show()