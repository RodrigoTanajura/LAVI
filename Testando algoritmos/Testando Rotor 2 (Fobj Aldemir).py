#%%
print("Iniciando bibliotecas e funções")
import time
st1 = time.time()
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
from scipy.optimize import differential_evolution, dual_annealing, shgo

st2 = time.time()
steel = rs.materials.steel
samples = 501
speed_range = np.linspace(0, 500, samples) # rads/s
kxx = 3e6
kyy = 6e6

cxx = 120
cyy = 100

bounds=[(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

# Rotor 2

shaft = [
        rs.ShaftElement(1.2 / (12), idl=0.03, odl=0.08, material=steel)
        for i in range(12)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(12 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
    ),
    rs.DiskElement.from_geometry(
        n=(7), material=steel, width=0.04, i_d=0.08, o_d=0.4
    ),
    rs.DiskElement.from_geometry(
        n=(5), material=steel, width=0.04, i_d=0.08, o_d=0.4
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
densidade_nos_picos = intervalo*2+1
picos = find_peaks(z1, prominence=5e-6)
n_de_picos = len(picos[0])
z2 = []
speed_range = []
for i in range(n_de_picos):
    speed_range = np.append(speed_range, speed[picos[0][i]-intervalo:picos[0][i]+intervalo])
    z2 = np.append(z2, z1[picos[0][i]-intervalo:picos[0][i]+intervalo])
plt.plot(speed,z1)
plt.show()
plt.plot(speed_range,z2,'o')

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
# Evolução diferencial
algoritmo = "Ev.Dif."
print("Começando otimização de evolução diferencial, tempo decorrido: {} segundos".format(round(et1-st1,2)))
for i in range(10):
    st = time.time()
    res = differential_evolution(objective, bounds, disp=True, polish=False)
    et = time.time()
    print(et-st)
    save_to_df()

et2 = time.time()
# print("Começando otimização de Dual Annealing, tempo decorrido: {} segundos".format(round(et2-et1)))
# # Dual annealing
# for i in range(10):
#     algoritmo = "D.Ann."
#     st = time.time()
#     res = dual_annealing(objective, bounds, maxfun=1000)
#     et = time.time()
#     save_to_df()

df1 = pd.DataFrame(df[0:10], columns=colunas)
# df2 = pd.DataFrame(df[10:20], columns=colunas)

with pd.ExcelWriter("Dados Rotor 2.xlsx", mode='w', engine='openpyxl') as writer:
    df1.to_excel(writer, sheet_name="Ev.Dif", index=False, float_format="%.2f")  
    # df2.to_excel(writer, sheet_name="D.Ann.", index=False, float_format="%.2f") 
# %%
