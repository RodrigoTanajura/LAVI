# lmfit
# <examples/doc_parameters_basic.py>
#%%
print("Iniciando bibliotecas e funções")
import time
st1 = time.time()
from matplotlib import pyplot as plt
import numpy as np
from random import randint, gauss
from math import trunc
from scipy.signal import find_peaks
import pandas as pd
import openpyxl
import ross as rs
import lmfit
from lmfit import Minimizer, Parameters, report_fit
# Create model
# Rotor 1
st2 = time.time()
steel = rs.materials.steel
samples = 501
speed_range = np.linspace(0, 500, samples) # rads/s

shaft = [
        rs.ShaftElement(1.2 / (6), idl=0.03, odl=0.08, material=steel)
        for i in range(6)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(6 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
    )
]

kxx = 8e5
kyy = 1e6

cxx = 12
cyy = 10

k = {'Kxx':kxx, 'Kyy':kyy, 'Cxx':cxx, 'Cyy':cyy}

bounds=[[trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)],
        [trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)],
        [trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)],
        [trunc(0.1*1.1*cyy),trunc(2*1.1*cyy)]]

x0 = np.array([(bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2])

params = Parameters()
params.add('Kxx', value=x0[0], min=bounds[0][0], max=bounds[0][1])
params.add('Kyy', value=x0[1], min=bounds[1][0], max=bounds[1][1])
params.add('Cxx', value=x0[2], min=bounds[2][0], max=bounds[2][1])
params.add('Cyy', value=x0[3], min=bounds[3][0], max=bounds[3][1])

def change_bearings(vals):
    try:
        vals = vals.valuesdict()
    except:
        k = [vals['Kxx'],vals['Kyy'],vals['Cyy'],vals['Cxx']]
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

# create data to be fitted
noise = []
gaussiano = 3e-7
for i in range(samples):
    noise.append(gauss(0,gaussiano))

# Gerando os dados pseudoexperimentais
def real_input(k):
    return change_bearings(k)+noise

z1 = real_input(k)

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


# define objective function: returns the array to be minimized
def objective(k):
    res = (z1 - change_bearings(k))
    return np.array(res)


# do fit, here with the default leastsq algorithm
minner = Minimizer(objective, params)
result = minner.minimize()

# write error report
report_fit(result)

# try to plot results
try:
    import matplotlib.pyplot as plt
    plt.plot(speed_range, z1, '+')
    plt.plot(speed_range, change_bearings(k))
    plt.show()
except ImportError:
    pass
#%%
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

algoritmo = "LM do lmfit"
print("Começando otimização de LM, tempo decorrido: {} segundos".format(round(et1-st1,2)))
n = 1
for i in range(n):
    st = time.time()
    res = least_squares(objective, x0, method='lm', verbose=1)
    et = time.time()
    save_to_df()

et2 = time.time()

df1 = pd.DataFrame(df[0:n], columns=colunas)

with pd.ExcelWriter("LM3.xlsx", mode='w', engine='openpyxl') as writer:
    df1.to_excel(writer, sheet_name="Least_squares", index=False, float_format="%.2f")  
print("Finalizado")

# %%
