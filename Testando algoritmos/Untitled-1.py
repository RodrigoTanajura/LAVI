#%%
# O objetivo é rodar uma série de otimizações e gravar todos os resultados. Diferentes algoritmos são rodados
# em loops.
print("Iniciando bibliotecas e funções")
import time
st1 = time.time()
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import norm 
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from math import trunc
from scipy.signal import find_peaks
import pandas as pd
import openpyxl
from scipy.optimize import differential_evolution, dual_annealing, shgo

st2 = time.time()
steel = rs.materials.steel
samples = 501 
speed_range = np.linspace(0, 500, samples) # rads/s
kxx = 8e5
kyy = 1e6
cxx = 12
cyy = 10

bounds=[(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

bounds1 = [(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy))]

bounds2 = [(trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

x0 = ((bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2 )

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

k = [kxx,kyy]
c = [cxx, cyy]

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

def change_stiffness(k):
    bearings = [
        rs.BearingElement(0, kxx=k[0], kyy=k[1],
        cxx=x0[2], cyy=x0[3]),
        rs.BearingElement(6, kxx=k[0], kyy=k[1],
        cxx=x0[2], cyy=x0[3])
    ]
    rotor = rs.Rotor(shaft, disks, bearings)
    results1 = rotor.run_freq_response(speed_range)
    r = results1.freq_resp[16][16]
    return np.sqrt(r.imag**2 + r.real**2)

def change_dampening(c):
    bearings = [
        rs.BearingElement(0, kxx=x1[0], kyy=x1[1],
        cxx=c[0], cyy=c[1]),
        rs.BearingElement(6, kxx=x1[0], kyy=x1[1],
        cxx=c[0], cyy=c[1]),
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
densidade_nos_picos = intervalo*2+1
picos = find_peaks(z1, prominence=5e-6)
n_de_picos = len(picos[0])
z2 = []
speed_range = []
for i in range(n_de_picos):
    speed_range = np.append(speed_range, speed[picos[0][i]-intervalo:picos[0][i]+intervalo])
    z2 = np.append(z2, z1[picos[0][i]-intervalo:picos[0][i]+intervalo])

z1 = z2

def objective_stiffness(k):
    res = (norm(z1 - change_stiffness(k), 2)**2)/(norm(z1, 2)**2)
    return res

def objective_dampening(k):
    res = (norm(z1 - change_dampening(c), 2)**2)/(norm(z1, 2)**2)
    return res

def objective(k):
    res = (norm(z1 - change_bearings(k), 2)**2)/(norm(z1, 2)**2)
    return res

# def objective_aldemir(k):
#     res = np.abs(change_bearings(k) - z1)
#     res = sum(res/np.abs(z1))
#     return res

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

# Otimizações:

res = differential_evolution(objective, bounds, disp=True, polish=False)

print('Tempo total:{}, Tempo 1:{}, Tempo 2:{}'.format(-st1+et2,et1-st1,et2-et1) )

print(res_dampening)
print(res)


# with pd.ExcelWriter("Dados Rotor 1-1.xlsx", mode='w', engine='openpyxl') as writer:
#     df1.to_excel(writer, sheet_name="Ev.Dif", index=False, float_format="%.2f")  
#     df2.to_excel(writer, sheet_name="D.Ann.", index=False, float_format="%.2f") 
# %%
