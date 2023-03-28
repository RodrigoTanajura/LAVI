#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import shgo
from math import trunc
from numpy import real, imag

#Testando SHGO no rotor do Aldemir.

freq = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 0]
amp = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 1]
im = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 2]

odl = 0.017 
z = np.sqrt(amp**2 + im**2)
z = z[25:,]
freq = freq[25:,]
speed_range = []
z1 = []
j = 0 
f_de_correcao = np.pi*2

reducao_de_densidade = 8
for i in range(int(len(freq)/reducao_de_densidade)):
    speed_range.append(freq[j])
    z1.append(z[j])
    j += reducao_de_densidade

speed_range = speed_range[0:50]
z1 = z1[0:50]

# k = [5e3,5e3,10,10,5e5,5e5,10,10]
k = [855130.682384417,52024576.8035698,7.45161127921376,25.5867100796101,1197651.61618692,702261906.673702, 33.6793402012332,91.0335488825413,2.73016153382567, 4.84982845688298e-06,770.441680165024]
material = rs.Material("Aldemir", 7850, E=205e9, G_s=None, Poisson=0.29, color='#525252')

L=[0, 25, 49, 65, 80, 110, 140, 161, 190, 220, 250, 285, 295, 305, 330, 360, 390, 
420, 450, 480, 510, 523, 533, 543, 570, 594, 630, 664, 700, 730, 760, 790, 830, 862]
L = [L[i] - L[i - 1] for i in range(1, len(L))]

list_odl = [odl]*len(L)
for i in [11,12,21,22]:
    list_odl[i] = list_odl[i]+0.020

list_odl[2] = 0.029
list_odl[3] = 0.029


disks = [rs.DiskElement(Id=0.007513248437500, Ip=0.003844540885417,
     m=2.6375, color='Firebrick', n=12),
    rs.DiskElement(Id=0.007547431937500, Ip=0.003862032635417,
     m=2.6495, color='Firebrick', n=22)
]

shaft = []
for i in range(len(L)):
    shaft.append(
        rs.ShaftElement(
            material=material,
            L=L[i]/1000,
            idl=0,
            odl=list_odl[i],
            rotary_inertia=False,
            shear_effects=False,
        ))

# Massas extras:

def width(massa):
    return (4*massa/(material.rho*np.pi*odl**2))

disks.append(rs.DiskElement.from_geometry(n=0, material=material, width=width(0.14945*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=1, material=material, width=width(0.14945*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=3, material=material, width=width(1.1252*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=30, material=material, width=width(0.0441*2), i_d=0, o_d=odl/10))

# Funcao que "remonta" o rotor.

def change_bearings(k):
    bearings = [
        rs.BearingElement(3, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3]),
        rs.BearingElement(30, kxx=k[4], kyy=k[5],
        cxx=k[6], cyy=k[7]),
    ]
    rotor = rs.Rotor(shaft, disks, bearings)
    results1 = rotor.run_freq_response(speed_range)
    results = results1.freq_resp[13*4][8*4]
    return np.sqrt(results.imag**2 + results.real**2)

def objective(k):
    res = np.abs(change_bearings(k) - z1)
    res = sum(res/np.abs(z1))
    return res

def checkin():
    print(f"Iteração ocorreu.")

res = shgo(objective, bounds=
[(5e5,5e6),(5e5,5e6),(0,2e2),(0,2e2),(5e5,5e6),(5e5,5e6),
(0,2e2),(0,2e2)],
callback=checkin())

# change_bearings(k)
print("Finalizado")
plt.plot(speed_range, z1,'b')
plt.plot(speed_range, change_bearings(k),'r')
plt.show()
# %%
