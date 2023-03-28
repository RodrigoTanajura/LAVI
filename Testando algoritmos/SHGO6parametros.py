#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import shgo
from math import trunc

steel = rs.materials.steel
samples = 50
speed_range = np.linspace(105, 135, samples) # rads/s
kx = 1e6*0.92435
ky = 1e6*1.14362
kxy = 1e6
cxx = 1e3*0.1234
cyy = 1e3*0.5748
cxy = 1e3

k = [kx,ky,kxy,cxx,cyy,cxy]

shaft = [
        rs.ShaftElement(1.2 / (6), idl=0.03, odl=0.08, material=steel)
        for i in range(6)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(6 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
    )
]


def change_bearings(k):
    bearings = [
        rs.BearingElement(0, kxx=k[0], kyy=k[1], kxy=k[2],
        cxx=k[3], cyy=k[4], cxy=k[5]),
        rs.BearingElement(6, kxx=k[0], kyy=k[1], kxy=k[2],
        cxx=k[3], cyy=k[4], cxy=k[5]),
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

z1 = real_input([kx,ky,kxy,cxx,cyy,cxy])

def objective(k):
    res = np.abs(change_bearings(k) - z1)
    res = sum(res/np.abs(z1))
    return res

res = shgo(objective, bounds=[(0.93*kx,1.17*kx),
(0.92*ky,1.18*ky),(0.95*kxy,1.12*kxy),(0.94*cxx,1.15*cxx),
(0.91*cyy,1.17*cyy),(0.91*cxy,1.18*cxy)], 
options = {"disp":"True","f_min":10})

# options= {"f_min":10e-6})

plt.plot(speed_range, change_bearings([kx,ky,kxy,cxx,cyy,cxy]),'b')
plt.plot(speed_range, change_bearings(res.x),'r')
plt.plot(speed_range, real_input([kx,ky,kxy,cxx,cyy,cxy]),'o')
plt.show()
print("K encontrado vs esperado:")
for i in range(len(res.x)):
    print (trunc(res.x[i]), "|" , trunc(k[i]))
print("Noise de sigma = ", gaussiano)

# %%
