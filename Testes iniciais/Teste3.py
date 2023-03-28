#%%
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from random import randint, gauss

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
t = 5
r = w/wn1
noise = []
noiseg = []
n = len(w)

for i in range(n):
    noise.append(gauss(0,0.4))

def xp(w, wn, csi):
    return (u/m)*((w/wn)**2)/((((1-(w/wn)**2)**2)+(2*csi*(w/wn))**2)**0.5)

parametros, pcov= curve_fit(xp, w,  xp(w, wn1, csi1)+noise, bounds = [0,np.inf])

errown = abs((wn1-parametros[0])/wn1)
errocsi = abs((csi1-parametros[1])/csi1)

plt.scatter(w,xp(w, wn1, csi1)+noise, label = 'FRF nominal')
y_line = xp(w,parametros[0],parametros[1])
plt.plot(w,y_line, label = 'FRF com noise gaussiano', color='r')
plt.legend(loc='best')
plt.show()
# print("Os parametros esperados eram wn = {0} e csi = {1}.".format(str(wn1), str(csi1)),'\n', "Os parametros obtidos foram wn* = {0} e csi* = {1}.".format(parametros[0], parametros[1]))
print("Erro em wn = (wn - wn*)/wn = {}%".format(round(errown,3)), '\n', "Erro em csi = (csi - csi*)/csi = {}%".format(round(errocsi, 3)))
print("Erro da matriz de Covariancia: wn = {}, csi = {}".format(pcov[0][0]**2,pcov[1][1]**2))
# %%
