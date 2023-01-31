#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import curve_fit, minimize_scalar, minimize, least_squares, shgo
from math import trunc

material = rs.Material("Aldemir", 7850, E=205e9, G_s=None, Poisson=0.29, color='#525252')
steel = rs.steel


samples = 101
speed_range = np.linspace(0, 500, samples) # rads/s
kx = 1e6*0.92435
ky = 1e6*1.14362
kxy = 1e6
cxx = 1e3*0.92435
cyy = 1e3*1.14362
cxy = 1e3

k = [kx,ky,kxy,cxx,cyy,cxy]

shaft = [
        rs.ShaftElement(0.86 / (33), idl=0, odl=0.017, material=material)
        for i in range(33)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(13), material=steel, width=0.020, i_d=0.017, o_d=0.150
    ),
    rs.DiskElement.from_geometry(
        n=(23), material=steel, width=0.020, i_d=0.017, o_d=0.150
    )
]
bearings = [
        rs.BearingElement(4, kxx=k[0], kyy=k[1], kxy=k[2],
        cxx=k[3], cyy=k[4], cxy=k[5]),
        rs.BearingElement(31, kxx=k[0], kyy=k[1], kxy=k[2],
        cxx=k[3], cyy=k[4], cxy=k[5])
]

rotor = rs.Rotor(shaft, disks, bearings)

# def change_bearings(k):
#     bearings = [
#         rs.BearingElement(4, kxx=k[0], kyy=k[1], kxy=k[2],
#         cxx=k[3], cyy=k[4], cxy=k[5]),
#         rs.BearingElement(31, kxx=k[0], kyy=k[1], kxy=k[2],
#         cxx=k[3], cyy=k[4], cxy=k[5]),
#     ]
#     rotor = rs.Rotor(shaft, disks, bearings)
#     results1 = rotor.run_freq_response(speed_range)
#     return np.abs(results1.freq_resp[16][16])

# noise = []
# gaussiano = 2e-6
# for i in range(samples):
#     noise.append(gauss(0,gaussiano))

# def real_input(k):
#     return change_bearings(k)

# def objective(k):
#     res = change_bearings(k) - real_input([kx,ky,kxy,cxx,cyy,cxy])
#     for i in range(len(res)):
#         res[i] = (res[i]**2)
#     res = sum(res)/2
#     return res

# res = shgo(objective, bounds=[(0.9*kx,1.1*kx),
# (0.9*ky,1.1*ky),(0.9*kxy,1.1*kxy),(0.9*cxx,1.1*cxx),(0.9*cyy,1.1*cyy),(0.9*cxy,1.1*cxy)])

# plt.plot(speed_range, change_bearings([kx,ky,kxy,cxx,cyy,cxy]),'b')
# plt.plot(speed_range, change_bearings(res.x),'r')
# plt.plot(speed_range, real_input([kx,ky,kxy,cxx,cyy,cxy]),'o')
# plt.show()
# print("K encontrado vs esperado:")
# for i in range(len(res.x)):
#     print (trunc(res.x[i]), "|", trunc(k[i]))
# print("Noise de sigma = ", gaussiano)
# %%
