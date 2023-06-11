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
from numpy.linalg import norm 
import openpyxl
from scipy.optimize import differential_evolution, dual_annealing, shgo, minimize
import lmfit
from lmfit import Minimizer, Parameters, report_fit
import scipy as sp

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

def FRF_only_stiffness(k):
    FRF_resp = np.zeros(len(speed_range))
    bearings = [
            rs.BearingElement(0, kxx=k[0], kyy=k[1],
            cxx=x0[2], cyy=x0[3]),
            rs.BearingElement(6, kxx=k[0], kyy=k[1],
            cxx=x0[2], cyy=x0[3])
        ]
    rotor = rs.Rotor(shaft, disks, bearings)
    Lambda, Phi = sp.linalg.eigh(rotor.K(0),rotor.M())
    Phi = Phi[:,:10]
    M_red = Phi.T@rotor.M()@Phi
    G_red = Phi.T@rotor.G()@Phi
    for i in range(len(speed_range)):
        w=speed_range[i]
        K_red = Phi.T @ rotor.K(w) @ Phi
        C_red = Phi.T @ rotor.C(w) @ Phi

        H_red = -w**2*M_red + (1j*w)*(C_red+w*G_red) + K_red
        Hinv_red = sp.linalg.inv(H_red)
        Hinv = Phi@Hinv_red@Phi.T
        FRF_resp[i] = np.abs(Hinv[16][16])
    return FRF_resp

def FRF_only_dampening(c):
    FRF_resp = np.zeros(len(speed_range))
    bearings = [
            rs.BearingElement(0, kxx=x1[0], kyy=x1[1],
            cxx=c[0], cyy=c[1]),
            rs.BearingElement(6, kxx=x1[0], kyy=x1[1],
            cxx=c[0], cyy=c[1])
        ]
    rotor = rs.Rotor(shaft, disks, bearings)
    Lambda, Phi = sp.linalg.eigh(rotor.K(0),rotor.M())
    Phi = Phi[:,:10]
    M_red = Phi.T@rotor.M()@Phi
    G_red = Phi.T@rotor.G()@Phi
    for i in range(len(speed_range)):
        w=speed_range[i]
        K_red = Phi.T @ rotor.K(w) @ Phi
        C_red = Phi.T @ rotor.C(w) @ Phi

        H_red = -w**2*M_red + (1j*w)*(C_red+w*G_red) + K_red
        Hinv_red = sp.linalg.inv(H_red)
        Hinv = Phi@Hinv_red@Phi.T
        FRF_resp[i] = np.abs(Hinv[16][16])
    return FRF_resp

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

# Funcao objetivo
def objective_stiffness(k):
    res = (norm(z1 - FRF_only_stiffness(k), 2)**2)/(norm(z1, 2)**2)
    return res

def objective_dampening(c):
    res = (norm(z1 - FRF_only_dampening(c), 2)**2)/(norm(z1, 2)**2)
    return res

def change_bearings(vals):
    try:
        vals = vals.valuesdict()
    except:
        pass
    k = [vals['Kxx'],vals['Kyy'],vals['Cyy'],vals['Cxx']]
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

def change_stiffness(vals):
    try:
        vals = vals.valuesdict()
    except:
        pass
    k = [vals['Kxx'],vals['Kyy']]
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

def change_dampening(vals):
    try:
        vals = vals.valuesdict()
    except:
        pass
    c = [vals['Cyy'],vals['Cxx']]
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

x0 = ((bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2 )

# Para salvar no Excel
df = []
              
colunas = ['Tempo','F_obj', 'kxx','Erro de kxx','kyy','Erro de kyy',
           'cxx','Erro de cxx','cyy','Erro de cyy','N_de_eval', "Algoritmo"]

def save_to_df():
    df.append((et-st,res_dampening.fun,trunc(res.x[0]),np.abs((k[0]-res.x[0])/(k[0]))*100,
    trunc(res.x[1]),np.abs((k[1]-res.x[1])/(k[1]))*100,
    trunc(res_dampening.x[0]),np.abs((k[0]-res_dampening.x[0])/(k[0]))*100,
    trunc(res_dampening.x[1]),np.abs((k[1]-res_dampening.x[1])/(k[1]))*100,
    res.nfev+res_dampening.fun,algoritmo))

n_de_testes = 1

# Otimizações:
# Evolução diferencial
algoritmo = "Ev.Dif."
print("Começando otimização de evolução diferencial.")
for i in range(n_de_testes):
    st = time.time()
    res = differential_evolution(objective_stiffness, bounds1, disp=False, polish=False)
    x1 = [res.x[0], res.x[1]]
    res_dampening = differential_evolution(objective_dampening, bounds2, disp=False, polish=True)
    et = time.time()
    save_to_df()

print("Começando otimização de Dual Annealing.")

# Dual annealing
for i in range(n_de_testes):
    algoritmo = "D.Ann."
    st = time.time()
    res = dual_annealing(objective_stiffness, bounds1, maxfun=1000)
    x1 = [res.x[0], res.x[1]]
    res_dampening = dual_annealing(objective_dampening, bounds2, maxfun=1000)
    et = time.time()
    save_to_df()

print("Começando otimização de Levenberg Marquardt.")

# Levenberg Marquardt
params1 = Parameters()
params1.add('Kxx', value=x0[0], min=bounds[0][0], max=bounds[0][1])
params1.add('Kyy', value=x0[1], min=bounds[1][0], max=bounds[1][1])

params2 = Parameters()
params2.add('Cxx', value=x0[2], min=bounds[2][0], max=bounds[2][1])
params2.add('Cyy', value=x0[3], min=bounds[3][0], max=bounds[3][1])

for i in range(n_de_testes):
    st = time.time()
    res = Minimizer(objective_stiffness, params1).minimize()
    x1 = [res.params.get('Kxx'), res.params.get('Kyy')]
    res_dampening = Minimizer(objective_dampening, params2).minimize()
    et = time.time()
    df.append((et-st, np.sqrt(sum(res_dampening.residual**2)), 
    trunc(res.params.get('Kxx')), np.abs((k[0]-res.params.get('Kxx'))/(k[0]))*100,
    trunc(res.params.get('Kxx')), np.abs((k[0]-res.params.get('Kxx'))/(k[0]))*100,
    trunc(res.params.get('Kxx')), np.abs((k[0]-res.params.get('Kxx'))/(k[0]))*100,
    trunc(res.params.get('Kxx')), np.abs((k[0]-res.params.get('Kxx'))/(k[0]))*100,
    res.nfev+res_dampening.nfev, algoritmo))

print("Começando otimização de Nelder Mead.")
# Nelder Mead
for i in range(n_de_testes):
    st = time.time()
    res = minimize(objective_stiffness, [x0[0],x0[1]], method=algoritmo, bounds=bounds1)
    x1 = [res.x[0], res.x[1]]
    print('Otimização do amortecimento')
    res_dampening = minimize(objective_dampening, [x0[2],x0[3]], method=algoritmo, bounds=bounds2)
    et = time.time()
    save_to_df()

df1 = pd.DataFrame(df[0:10], columns=colunas)
df2 = pd.DataFrame(df[10:20], columns=colunas)
df3 = pd.DataFrame(df[20:30], columns=colunas)
df4 = pd.DataFrame(df[30:40], columns=colunas)

with pd.ExcelWriter("Dados Rotor 1.xlsx", mode='w', engine='openpyxl') as writer:
    df1.to_excel(writer, sheet_name="Ev.Dif", index=False, float_format="%.2f")  
    df2.to_excel(writer, sheet_name="D.Ann.", index=False, float_format="%.2f")
    df3.to_excel(writer, sheet_name="LM", index=False, float_format="%.2f") 
    df4.to_excel(writer, sheet_name="NM", index=False, float_format="%.2f")  
# %%
