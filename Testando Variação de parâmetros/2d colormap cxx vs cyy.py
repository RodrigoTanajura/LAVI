# GRID 10x10 para kxx vs kyy
#%%
import ross as rs
import numpy as np
from numpy.linalg import norm
from scipy.signal import find_peaks
from random import randint, gauss
from math import trunc
from matplotlib import pyplot as plt
from matplotlib import cm

# Definindo o rotor e suas características
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

def change_bearings(cx, cy):
    bearings = [
        rs.BearingElement(0, kxx=kxx, kyy=kyy,
        cxx=cx, cyy=cy),
        rs.BearingElement(6, kxx=kxx, kyy=kyy,
        cxx=cx, cyy=cy)
    ]
    rotor = rs.Rotor(shaft, disks, bearings)
    results1 = rotor.run_freq_response(speed_range)
    r = results1.freq_resp[16][16]
    return np.sqrt(r.imag**2 + r.real**2)

noise = []
gaussiano = 0
for i in range(samples):
    noise.append(gauss(0,gaussiano))

# Gerando os dados pseudoexperimentais
def real_input(kx, ky):
    return change_bearings(kx, ky)+noise

z1 = real_input(cxx, cyy)

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

# Definindo a função objetivo.
def objective(cx, cy):
    res = (norm(z1 - change_bearings(cx, cy), 2)**2)/(norm(z1, 2)**2)
    return res

# n é o número de valores de kx e ky. n = 10 irá gerar um grid de 100 pontos, n = 20 irá gerar um grid de 400 pontos e assim por diante.
# Cada ponto é uma avaliação da função objetivo.
n = 20   #obs: n = 10 demora na ordem de 60 segundos para executar.

# 0.6545 é o tempo de realizar uma avaliação da função objetivo (no meu computador).
print("Tempo estimado:{}".format((n**2)*0.6545))

# Gerando os valores de x e y respectivamente que serão inputados.
x = np.linspace(bounds[2][0], bounds[2][1], n)
y = np.linspace(bounds[3][0], bounds[3][1], n)

# l será a lista de pares de kx e ky que serão inputados na função objetivo
l = []
for i in x:
    for j in y:
        l.append([i,j])

# f(x,y) é um array com as avaliações de todos os pares pela função objetivo.
def f(x, y):
    return np.array([objective(l[i][0],l[i][1]) for i in range(len(l))])

# Armazenando os valores das avaliações da função objetivo em formato de matriz n,n.
z = f(x, y).reshape(n,n)

# A matrix precisa ser transposta senão kx e ky acabam trocados nos eixos x e y do gráfico. 
z = np.matrix.transpose(z)

c = plt.imshow(z, interpolation='bilinear', extent=[bounds[2][0], bounds[2][1], bounds[3][0], bounds[3][1],])
colorbar = plt.colorbar(c)
colorbar.set_label('Função objetivo')
plt.xlabel("cxx")
plt.ylabel("cyy")
plt.show()
# %%
