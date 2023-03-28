#%%
# NSGA 2
print("Iniciando bibliotecas e funções")
import time
st1 = time.time()
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from math import trunc
from scipy.signal import find_peaks
import pandas as pd
import openpyxl

st2 = time.time()
steel = rs.materials.steel
samples = 501
speed_range = np.linspace(0, 500, samples) # rads/s
kxx = 8e5
kyy = 1e6

cxx = 12
cyy = 10

bounds=[(0.8*1.1*kxx,1.2*1.1*kxx),
        (0.8*1.1*kyy,1.2*1.1*kyy),
        (0.1*1.1*cxx,2*1.1*cxx),
        (0,1*1.1*cyy,2*1.1*cyy)]

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

speed = speed_range

# Localizando picos e alterando o speed range pra focar neles.
intervalo = 10
densidade_nos_picos = intervalo*2
picos = find_peaks(z1, prominence=5e-6)
n_de_picos = len(picos[0])
z2 = []
speed_range = []
for i in range(n_de_picos):
    speed_range = np.append(speed_range, np.linspace(picos[0][i]-intervalo,picos[0][i]+intervalo,densidade_nos_picos))
    z2 = np.append(z2, z1[picos[0][i]-intervalo:picos[0][i]+intervalo])

z1 = z2

# Funcao objetivo
def objective(k):
    res = np.abs(change_bearings(k) - z1)
    res = sum(res/np.abs(z1))
    return res

x0 = ((bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2 )

# Para salvar no Excel
df = []
              
colunas = ['Tempo','F_obj', 'kxx','Erro de kxx','kyy','Erro de kyy',
           'cxx','Erro de cxx','cyy','Erro de cyy','N_de_eval', "Algoritmo"]

def save_to_df():
    df.append((et-st,res.fun,trunc(res.x[0]),np.abs((k[0]-res.x[0])/(k[0]))*100,
    trunc(res.x[1]),np.abs((k[1]-res.x[1])/(k[1]))*100,
    trunc(res.x[2]),np.abs((k[2]-res.x[2])/(k[2]))*100,
    trunc(res.x[3]),np.abs((k[3]-res.x[3])/(k[3]))*100,
    res.nfev,algoritmo))
et1 = time.time()

# Otimizações:
algoritmo = "NSGA2"
algorithm = NSGA2(pop_size=100)

res = minimize(objective,
               algorithm,
               ('n_gen', 200),
               verbose=True)

# %%
