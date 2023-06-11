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

def change_bearings(kx, ky):
    bearings = [
        rs.BearingElement(0, kxx=kx, kyy=ky,
        cxx=cxx, cyy=cyy),
        rs.BearingElement(6, kxx=kx, kyy=ky,
        cxx=cxx, cyy=cyy)
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

z1 = real_input(kxx, kyy)

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
def objective(kx, ky):
    res = (norm(z1 - change_bearings(kx, ky), 2)**2)/(norm(z1, 2)**2)
    return res

n = 10 

x = np.linspace(bounds[0][0], bounds[0][1], n)
y = np.linspace(bounds[1][0], bounds[1][1], n)

g = np.array([objective(i,kyy) for i in x])
c = np.array([objective(kxx, j) for j in y])

plt.plot(np.arange(bounds[0][0],bounds[0][1], 352*100), g)
plt.title("Variando somente kxx")
plt.xlabel("Kxx variando de {} a {}".format(bounds[0][0],bounds[0][1]))
plt.ylabel("Função objetivo")
plt.show()
plt.plot(np.arange(bounds[1][0],bounds[1][1], 440*100), c)
plt.title("Variando somente kyy")
plt.xlabel("Kyy variando de {} a {}".format(bounds[1][0],bounds[1][1]))
plt.ylabel("Função objetivo")
plt.show()
# %%
