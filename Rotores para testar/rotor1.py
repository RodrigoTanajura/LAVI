#%%
import ross as rs
from math import trunc
import numpy as np

steel = rs.materials.steel

kxx = 8e5
kyy = 1e6

cxx = 12
cyy = 10

k = [kxx,kyy,cxx,cyy]

bounds=[(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

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

rotor1 = rs.Rotor(shaft, disks, bearings)
rotor1.plot_rotor()

# %%
