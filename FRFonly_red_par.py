# %%
"""
Created on Wed Apr 12 21:07:07 2023

@author: Guilherme
"""
import multiprocessing as mp
import ipyparallel as ipp
import pandas as pd
import time
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import norm 
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from math import trunc
from scipy.signal import find_peaks
import openpyxl
import scipy as sp
from scipy.optimize import minimize, least_squares
from scipy.optimize import differential_evolution, dual_annealing, shgo


rc = ipp.Cluster(n=mp.cpu_count()).start_and_connect_sync()
dv = rc[:]
#rc.ids

dv.block = True
dv.activate('')

#speed_range = np.linspace(0, 2500, samples)

steel = rs.materials.steel
samples = 501 
speed_range = np.linspace(0, 500, samples) # rads/s
kxx = 1e6
kyy = 1e6
cxx = 18
cyy = 5

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

noise = []
gaussiano = 8e-7
for i in range(samples):
    noise.append(gauss(0,gaussiano))

# Gerando os dados pseudoexperimentais

def real_input(k):
    return change_bearings(k)+noise

z1 = real_input([kxx,kyy,cxx,cyy]) # z1 é o vetor de amplitudes dos dados experimentais, de 0 a 500.

speed = speed_range 

intervalo = 10
densidade_nos_picos = intervalo*2+1
picos = find_peaks(z1, prominence=5e-6)
n_de_picos = len(picos[0])
z2 = []
speed_range = []
for i in range(n_de_picos):
    speed_range = np.append(speed_range, speed[picos[0][i]-intervalo:picos[0][i]+intervalo])
    z2 = np.append(z2, z1[picos[0][i]-intervalo:picos[0][i]+intervalo])

bounds1 = [(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy))]

z1 = z2 
# z1 agora é as amplitudes somente no entorno dos picos.
# speed_range é agora somente 60 pontos.

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

rc = ipp.Cluster(n=mp.cpu_count()).start_and_connect_sync()
dv = rc[:]
#rc.ids

dv.execute("""
""")

dv.block = True
dv.activate('')

dv.execute("""
steel = rs.materials.steel
samples = 501 
speed_range = np.linspace(0, 500, samples) # rads/s
kxx = 1e6
kyy = 1e6
cxx = 18
cyy = 5

shaft = [
    rs.ShaftElement(1.2 / (6), idl=0.03, odl=0.08, material=steel)
    for i in range(6)
]

disks = [
rs.DiskElement.from_geometry(
    n=(6 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
)
]

x0 = ((bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2 )
""")

def parallelization(k):
    dv.apply_sync(FRF_only_stiffness, k)
    dv.scatter("z1", z1, broadcast=True)
    dv.scatter("speed_range",speed_range)
    dv.execute("""
    df = FRF_only_stiffness(k)
    """)
    return dv.gather(["df"])

def objective(k):
    return parallelization(k)

print(objective(k))

# res = differential_evolution(objective, bounds1, disp=True, polish=False)


# %%
