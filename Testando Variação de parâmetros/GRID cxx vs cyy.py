# GRID 20x20 para cxx vs cyy
#%%
import ross as rs
import numpy as np
from numpy.linalg import norm
from scipy.signal import find_peaks
from random import randint, gauss
from math import trunc
from matplotlib import pyplot as plt
from matplotlib import cm

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
gaussiano = 3e-7
for i in range(samples):
    noise.append(gauss(0,gaussiano))

# Gerando os dados pseudoexperimentais
def real_input(cx, cy):
    return change_bearings(cx, cy)+noise

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

def objective(cx, cy):
    res = (norm(z1 - change_bearings(cx, cy), 2)**2)/(norm(z1, 2)**2)
    return res

n = 20
print("Tempo estimado:{}".format((n**2)*0.6545))

x = np.linspace(bounds[2][0], bounds[2][1], n)
y = np.linspace(bounds[3][0], bounds[3][1], n)

X, Y = np.meshgrid(x, y)

Z =[]

for j in range(n):
    Z.append([])
    for i in range(n):
        Z[j].append(objective(x[i],y[j]))
Z = np.array(Z)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, Y, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_xlabel('cxx')
ax.set_ylabel('cyy')
ax.set_zlabel('F_obj')
# %%
