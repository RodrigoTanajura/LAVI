#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import differential_evolution
from math import trunc
from numpy import real, imag
import time

steel = rs.materials.steel
samples = 50
speed_range = np.linspace(105, 135, samples) # rads/s
kxx = 1e6*0.92435
kyy = 1e6*1.14362
kxy = 1e6
cxx = 1e3*0.1234
cyy = 1e3*0.5748
cxy = 1e3

k = [kxx,kyy,cxx,cyy]

shaft = [
        rs.ShaftElement(1.2 / (6), idl=0.03, odl=0.08, material=steel)
        for i in range(6)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(6 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
    )
]

rotor = rs.Rotor(shaft, disks)

def change_bearings(k):
    bearings = [
        rs.BearingElement(0, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3]),
        rs.BearingElement(6, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3])
    ]
    rotor = rs.Rotor(shaft, disks, bearings)
    results1 = rotor.run_freq_response(speed_range)
    r = results1.freq_resp[16][16]
    return np.sqrt(r.imag**2 + r.real**2)

noise = []
gaussiano = 1e-6
for i in range(samples):
    noise.append(gauss(0,gaussiano))

def real_input(k):
    return change_bearings(k)+noise

z1 = real_input([kxx,kyy,cxx,cyy])

def objective(k):
    res = np.abs(change_bearings(k) - z1)
    res = sum(res/np.abs(z1))
    return res

bounds=[(0.93*kxx,1.17*kxx),
(0.92*kyy,1.18*kyy),(0.94*cxx,1.15*cxx),
(0.91*cyy,1.17*cyy)]

x0 = ((0.93*kxx+1.17*kxx)/2,(0.92*kyy+1.18*kyy)/2,
(0.94*cxx+1.15*cxx)/2, (0.91*cyy+1.17*cyy)/2 )

st = time.time()
res = differential_evolution(objective, bounds, disp=True)
et = time.time()

plt.plot(speed_range, change_bearings([kxx,kyy,cxx,cyy]),'b')
plt.plot(speed_range, change_bearings(res.x),'r')
plt.plot(speed_range, real_input([kxx,kyy,cxx,cyy]))
plt.show()
print("K encontrado | esperado | erro(%):")
for i in range(len(res.x)):
    print (trunc(res.x[i]), "|" , trunc(k[i]),"|",
            np.abs((k[i]-res.x[i])/(k[i])*100))
print("Noise de sigma = ", gaussiano,'\n', 'Tempo de calculo',trunc(et-st))
# %%
