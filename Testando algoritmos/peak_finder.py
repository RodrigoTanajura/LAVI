#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from math import trunc
from scipy.signal import find_peaks
import pandas as pd
import time
import openpyxl
from scipy.optimize import differential_evolution, dual_annealing, shgo

steel = rs.materials.steel
samples = 501
speed_range = np.linspace(0, 500, samples) # rads/s
kxx = 8e5
kyy = 1e6

cxx = 12
cyy = 10
n=6
# Rotor 1

shaft = [
        rs.ShaftElement(1.2 / (6), idl=0.03, odl=0.08, material=steel)
        for i in range(6)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(6 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
    )
]

k = [kxx,kyy,cxx,cyy]

bearings = [
        rs.BearingElement(0, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3]),
        rs.BearingElement(6, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3])
    ]

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
gaussiano = 3e-7
for i in range(samples):
    noise.append(gauss(0,gaussiano))

# Gerando os dados pseudoexperimentais
def real_input(k):
    return change_bearings(k)+noise

z1 = real_input([kxx,kyy,cxx,cyy])

# plt.plot(speed_range,z1)
# # Localizando picos
# st = time.time()
# densidade_nos_picos = 50
# picos = find_peaks(z1, prominence=5e-6)
# speed_range = np.linspace(picos[0][0]-5,picos[0][0]+5,densidade_nos_picos)
# speed_range = np.append(speed_range, np.linspace(picos[0][1]-10,picos[0][1]+10,densidade_nos_picos))
# et = time.time()
# print(et-st)
# plt.plot(speed_range, change_bearings(k),'o')
# versao antiga acima ^
exp_data = z1
peak_prominence = 5e-6
peak_interval = 10

peaks = find_peaks(exp_data, prominence=peak_prominence)
peak_speed_range = []
peak_exp_data = []
for i in range(len(peaks[0])):
    peak_speed_range = np.append(peak_speed_range, speed_range[peaks[0][i]-peak_interval:peaks[0][i]+peak_interval])
    peak_exp_data = np.append(peak_exp_data, exp_data[peaks[0][i]-peak_interval:peaks[0][i]+peak_interval])

plt.plot(speed_range, z1)
plt.show()

plt.plot(peak_speed_range,peak_exp_data, 'o')
plt.show()
# %%