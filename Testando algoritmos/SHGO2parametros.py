#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import curve_fit, minimize_scalar, minimize, least_squares, shgo
   

steel = rs.materials.steel
samples = 101
speed_range = np.linspace(0, 500, samples) # rads/s
kx = 1e6*0.924354642354235
ky = 1e6*1.1436286534224323

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
        rs.BearingElement(0, kxx=k[0], kyy=k[1], cxx=0),
        rs.BearingElement(6, kxx=k[0], kyy=k[1], cxx=0),
    ]
    rotor = rs.Rotor(shaft, disks, bearings)
    results1 = rotor.run_freq_response(speed_range)
    return np.abs(results1.freq_resp[16][16])

noise = []
for i in range(samples):
    noise.append(gauss(0,5e-6))

def real_input(k):
    return change_bearings(k) + noise

def objective(k):
    res = change_bearings(k) - real_input([kx,ky])
    for i in range(len(res)):
        res[i] = (res[i]**2)
    res = sum(res)/2
    return res

# res = least_squares(objective, x0=8e5, verbose=2,
#  method='trf', xtol=1e-8, ftol=1e-8, gtol = 1e-3)

res = shgo(objective, bounds=[(0.9*kx,1.1*kx),
(0.9*ky,1.1*ky)])

plt.plot(speed_range, change_bearings([kx, ky]),'b')
plt.plot(speed_range, change_bearings(res.x),'r')
plt.plot(speed_range, real_input([kx, ky]),'o')
plt.show()
print("K encontrado foi:",res.x[0], "e", res.x[1])
res.x[0]==kx
res.x[1]==ky
# %%
