# Version 2 of parameter identification tool. Dividing the process 
# into two steps: stiffness identification and dampening identification.
#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
from numpy.linalg import norm 
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import differential_evolution, dual_annealing
from math import trunc
from numpy import real, imag
import pandas as pd
import time
import openpyxl
from scipy.signal import find_peaks

steel = rs.materials.steel
samples = 501
speed_range = np.linspace(0, 500, samples) # rads/s
gaussiano = 3e-7
kxx = 1e6*0.92435
kyy = 1e6*1.14362
kxy = 1e6
cxx = 1e3*0.1234
cyy = 1e3*0.5748
cxy = 1e3

k = [kxx,kyy,cxx,cyy,kxx,kyy,cxx,cyy]

shaft = [
        rs.ShaftElement(1.2 / (6), idl=0.03, odl=0.08, material=steel)
        for i in range(6)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(6 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
    )
]
bearings = [
    rs.BearingElement(0, kxx=k[0], kyy=k[1],
    cxx=k[2], cyy=k[3]),
    rs.BearingElement(6, kxx=k[0], kyy=k[1],
    cxx=k[2], cyy=k[3])
]

rotor = rs.Rotor(shaft, disks, bearings)

results1 = rotor.run_freq_response(speed_range)
r = results1.freq_resp[16][16]

def noise():
    noise = []
    for i in range(samples):
        noise.append(gauss(0,gaussiano))
    return noise

exp_data = np.sqrt(r.imag**2 + r.real**2) + noise()

bounds=[(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy)),
        (trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

x0 = ((bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2 )

b_pos = [0,6]

# Acrescentar: kwargs para alterar parametros de otimização, possibilidade de rodar a simulação várias vezes e tomar a média
# opção de acionar valores diagonais nas matrizes de rigidez e amortecimento
# 
def find_parameters(exp_data, speed_range, bounds, b_pos, b_p, inp, out, method ='LM',  objective=1, diag_values=False, peak_interval=10, peak_prominence=5e-6):
    # Finding peaks
    peaks = find_peaks(exp_data, prominence=peak_prominence)
    peak_speed_range = []
    peak_exp_data = []
    for i in range(len(peaks[0])):
        peak_speed_range = np.append(peak_speed_range, speed_range[peaks[0][i]-peak_interval:peaks[0][i]+peak_interval])
        peak_exp_data = np.append(peak_exp_data, exp_data[peaks[0][i]-peak_interval:peaks[0][i]+peak_interval])
    # Changing bearings and calling the rotor function.
    def change_bearings(k):
        bearings = []
        for i in range(int(len(b_p)/4)):
            bearings.append(rs.BearingElement(b_pos[i], kxx=k[0+(i*4)], kyy=k[1+(i*4)],
                                            cxx=k[2+(i*4)], cyy=k[3+(i*4)])),
        return rs.Rotor(shaft, disks, bearings)
    def tentativa(k):
        change_bearings(k)
        results1 = rotor.run_freq_response(peak_speed_range)
        r = results1.freq_resp[inp][out]
        return np.sqrt(r.imag**2 + r.real**2)
    if objective == 1 and method != 'LM':
        def obj(k):
            res = (norm(peak_exp_data - tentativa(k), 2)**2)/(norm(peak_exp_data, 2)**2)
            return res
    elif objective == 2 and method != 'LM':
        def obj(k):
            res = np.abs(tentativa(k) - peak_exp_data)
            res = sum(res/np.abs(peak_exp_data))
            return res
    if method == 'LM':
        print('LM')
    elif method == 'diff_ev':
        st = time.time()
        res = differential_evolution(obj, bounds, disp=True)
        et = time.time()
        return res.x, et-st
    
# res = find_parameters(exp_data, speed_range, bounds, b_pos, k, 4, 4, method='diff_ev')

# print ("Optimized bearing parameters\n",
#                "bearing_number | kxx | kyy | cxx | cyy\n")
# for i in range(2):
#     print("{0} | {1} | {2} | {3} | {4}\n".format(i,res.x[0+(i*4)],res.x[1+(i*4)],res.x[2+(i*4)],res.x[3+(i*4)]))    
# print("\n Run time:{}\n".format(et-st),"Number of evaluations:{}".format(res.evals))

# colunas = ['Tempo','F_obj', 'kxx','Erro de kxx','kyy','Erro de kyy','cxx','Erro de cxx','cyy','Erro de cyy','N_de_eval']

# def save_to_df():
#     df.append((et-st,res.fun,trunc(res.x[0]),np.abs((k[0]-res.x[0])/(k[0]))*100,
#     trunc(res.x[1]),np.abs((k[1]-res.x[1])/(k[1]))*100,
#     trunc(res.x[2]),np.abs((k[2]-res.x[2])/(k[2]))*100,
#     trunc(res.x[3]),np.abs((k[3]-res.x[3])/(k[3]))*100,
#     res.nfev))

# for i in range(10):
#     st = time.time()
#     res = differential_evolution(objective, bounds, disp=True)
#     et = time.time()
#     save_to_df()

# data = pd.DataFrame(df, columns=colunas)
# data.to_excel("Testando Algoritmos.xlsx", index=False)  

# %%
