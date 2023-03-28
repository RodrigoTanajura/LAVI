# Testando algoritmo de dual annealing
#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import differential_evolution, dual_annealing
from math import trunc
from numpy import real, imag
import pandas as pd
import time
import openpyxl

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

df = []
              
colunas = ['Tempo','F_obj', 'kxx','Erro de kxx','kyy','Erro de kyy','cxx','Erro de cxx','cyy','Erro de cyy','N_de_eval']

def save_to_df():
    df.append((et-st,res.fun,trunc(res.x[0]),np.abs((k[0]-res.x[0])/(k[0]))*100,
    trunc(res.x[1]),np.abs((k[1]-res.x[1])/(k[1]))*100,
    trunc(res.x[2]),np.abs((k[2]-res.x[2])/(k[2]))*100,
    trunc(res.x[3]),np.abs((k[3]-res.x[3])/(k[3]))*100,
    res.nfev))

for i in range(3):
    st = time.time()
    res = dual_annealing(objective, bounds, maxfun=500)
    et = time.time()
    save_to_df()

data = pd.DataFrame(df, columns=colunas)
data.to_excel("Testando Algoritmos 2.xlsx", index=False) 
# %%
