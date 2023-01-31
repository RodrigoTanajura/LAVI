#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import shgo
from math import trunc

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

k = [8.551e5,1.198e6,7.452,33.679,5.202e7,7.023e8,25.587,91.033,2.730,4.85e-6]
material = rs.Material("Aldemir", 7850, E=205e9, G_s=None, Poisson=0.29)
steel = rs.materials.steel
odl = 0.017 
samples = 97
speed_range = freq # rads/s

disks = [
    rs.DiskElement6DoF.from_geometry(
        n=(13), material=steel, width=0.020, i_d=0.017, o_d=0.150
    ),
    rs.DiskElement6DoF.from_geometry(
        n=(23), material=steel, width=0.020, i_d=0.017, o_d=0.150
    )
]

bearings = [
        rs.BearingElement6DoF(4, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3]),
        rs.BearingElement6DoF(31, kxx=k[4], kyy=k[5],
        cxx=k[6], cyy=k[7]),
    ]
    
shaft = [
    rs.ShaftElement6DoF(
        material=material,
        L=(0.86/33),
        idl=0,
        odl=odl,
        alpha= k[8],
        beta=k[9],
        rotary_inertia=False,
        shear_effects=False,
    )
    for i in range(33)
]

rotor = rs.Rotor(shaft, disks, bearings)
results1 = rotor.run_freq_response(speed_range)
results1.plot_magnitude(inp = 48, out = 48)
# %%
