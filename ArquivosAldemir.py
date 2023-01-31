#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import shgo
from math import trunc

# Opening the files to define de frequency range vector as well as the amplitudes vector
f = open(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT", "r")

freq = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 0]
amp = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 1]
im = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 2]
z = []
f_de_correcao = np.pi*0.017
freq = freq*f_de_correcao

for i in range(len(freq)):
    z.append(np.abs((amp[i]**2 + im[i]**2)**(0.5)))

f.close()


# k = [5e5,5e5,1e2,1e2,5e7,5e7,1e2,1e2]
k = [8.551e5,1.198e6,7.452,33.679,5.202e7,7.023e8,25.587,91.033]
material = rs.Material("Aldemir", 7850, E=205e9, G_s=None, Poisson=0.29, color='#525252')
steel = rs.materials.steel
odl = 0.017
samples = 97
speed_range = freq # rads/s

disks = [
    rs.DiskElement.from_geometry(
        n=(13), material=steel, width=0.020, i_d=0.017, o_d=0.150
    ),
    rs.DiskElement.from_geometry(
        n=(23), material=steel, width=0.020, i_d=0.017, o_d=0.150
    )
]

def change_bearings(k):
    bearings = [
        rs.BearingElement(4, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3]),
        rs.BearingElement(31, kxx=k[4], kyy=k[5],
        cxx=k[6], cyy=k[7]),
    ]
    shaft = [
    rs.ShaftElement(
        material=material,
        L=(0.86/33),
        idl=0,
        odl=odl,
        rotary_inertia=True,
        shear_effects=True,
    )
    for i in range(33)
    ]
    rotor = rs.Rotor(shaft, disks, bearings)
    results1 = rotor.run_freq_response(speed_range)
    return np.abs(results1.freq_resp[48][48])

def objective(k):
    res = change_bearings(k) - z
    for i in range(len(res)):
        res[i] = (res[i]**2)
    res = sum(res)/2
    return res

def checkin():
    print(f"Iteração ocorreu.")

# res = shgo(objective, bounds=
# [(5e5,5e6),(5e5,5e7),(0,2e2),(0,2e2),(5e7,1e9),(5e7,1e9),(0,2e2),(0,2e2),(0,5),(0,1e-5)],
# callback=checkin())

# plt.plot(speed_range, z,'b')
plt.plot(speed_range, change_bearings(k),'r')
plt.show()


# %%
