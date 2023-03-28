#%%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

#parametros
x0 = 1
v0 = 1
w = np.linspace(0,100,100)
k = 10000
c = 100
m = 10
wn1 = ((k/m)**(1/2))
u = 10
csi1 = c/(2*m*wn1)
t = 10
r = w/wn1

def xp(w, wn, csi):
    return (u/m)*((w/wn)**2)/((((1-(w/wn)**2)**2)+(2*csi*(w/wn))**2)**0.5)+np.randon(0,0.1)

parametros, pcov= curve_fit(xp, w,  xp(w, wn1, csi1))

plt.scatter(w,xp(w, wn1, csi1), label = 'FRF nominal')
y_line = xp(w,parametros[0],parametros[1])
plt.plot(w,y_line, label = 'FRF ajustada', color='r')
plt.legend(loc='best')
plt.show()
print('Os parametros esperados eram', wn1, csi1, '\n', 'Os parametros encontrados foram:', parametros)

# %%
