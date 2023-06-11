#%%
# Testando otimizacao em duas etapas. Primeiro buscamos encontrar as frequencias naturais do sistema (rigidezes) e 
# depois refazemos a otimizacao para ajustar a amplitude dos picos (amortecimento).
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
import scipy as sp
from scipy.optimize import minimize, least_squares
from scipy.optimize import differential_evolution, dual_annealing, shgo

#Definindo o rotor de teste. Rotor simples de 6 nós com 1 disco central.

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

# Definindo quais sao os bounds do problema. Ou seja, qual a área de busca de cada variável.
# bounds são as 4 variáveis? kxx, kyy, cxx, cyy
# bounds1 são somente as rigidezes.
# bounds2 são somente os amortecimentos.

bounds=[(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

bounds1 = [(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy))]

bounds2 = [(trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

# X0 é o vetor que marca o ponto médio das condicoes de contorno.

x0 = ((bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2 )

# As funcoes change servem o proposito de redefinir os coeficientes dos mancais.
# Ou seja, sao as funcoes que serao chamadas pela funcao objetivo para fazer novas tentativas de parametros.
# Ha tres funcoes change. Se quisermos alterar todos os parametros simultaneamente usamos change_bearings.
# Se quisermos alterar somente as rigidezes ou somente os amortecimentos podemos usar as outras duas seguintes.
# Isto sera util para otimizar os parametros em etapas separadas.

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

def FRF_only_stiffness(k, pos):
    FRF_resp = np.zeros(len(speed_range))
    bearings = [
            rs.BearingElement(pos[i], kxx=k[0], kyy=k[1],
            cxx=x0[2], cyy=x0[3]),
            rs.BearingElement(pos[i], kxx=k[0], kyy=k[1],
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
# Noise é o ruído acrescentado aos dados gerados pelo ross para simular as incertezas experimentais.

noise = []
gaussiano = 8e-7
for i in range(samples):
    noise.append(gauss(0,gaussiano))

# Gerando os dados pseudoexperimentais

def real_input(k):
    return change_bearings(k)+noise

z1 = real_input([kxx,kyy,cxx,cyy]) # z1 é o vetor de amplitudes dos dados experimentais, de 0 a 500.

speed = speed_range #speed é o vetor de frequencias

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

z1 = z2 # z1 agora é as amplitudes somente no entorno dos picos.
z = []
# Funcao objetivo Ritto
def objective_stiffness(k, pos):
    res = []
    for i in range(len(z)):
        res.append((norm(z[i] - FRF_only_stiffness(k, pos[i]), 2)**2)/(norm(z[i], 2)**2))
    return sum(res)

def objective_dampening(c):
    res = (norm(z1 - FRF_only_dampening(c), 2)**2)/(norm(z1, 2)**2)
    return res
# return da funcao objetivo para amortecimento é multiplicado por 100, para aumentar a sensibilidade da resposta. Ou seja, a funcao objetivo terá de ser minimizada 100 vezes mais que para a rigidez.

# Otimizações:
# Somente as rigidezes:
st1 = time.time()
res = differential_evolution(objective_stiffness, bounds1, disp=True, polish=False)
et1 = time.time()

x1 = [res.x[0], res.x[1]] # Fixando as rigidezes com o resultado da otimizacao acima. Esses dados serao usados a seguir na otimizacao do amortecimento. 

print("Otimizando amortecimento")
# Agora somente o amortecimento:
res_dampening = differential_evolution(objective_dampening, bounds2, disp=True, polish=True)
et2 = time.time()

print('Tempo total:{}\n Tempo 1:{}\n Tempo 2:{}\n'.format(-st1+et2,et1-st1,et2-et1) )

print("Valores de rigidez:{}\n".format(res.x))
print("Valores de amortecimento:{}".format(res_dampening.x))

# %%
