#%%
import ross as rs
import numpy as np
from scipy import signal

steel = rs.materials.steel
kxx = 1.25e6
kyy = 1.4e6

cxx = 1e2
cyy = 180


k = [kxx,kyy,cxx,cyy]

shaft = [
        rs.ShaftElement(1.2 / (6), idl=0.03, odl=0.08, material=steel)
        for i in range(6)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(6 / 2), material=steel, width=0.08, i_d=0.08, o_d=0.4
    ),
    rs.DiskElement.from_geometry(
        n=(4), material=steel, width=0.04, i_d=0.08, o_d=0.4
    ),
    rs.DiskElement.from_geometry(
        n=(2), material=steel, width=0.04, i_d=0.08, o_d=0.4
    )
]

bearings = [
        rs.BearingElement(0, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3]),
        rs.BearingElement(6, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3])
    ]

speed_range = np.linspace(0,500,501)

rotor2 = rs.Rotor(shaft, disks, bearings)
rotor2.plot_rotor()
# results1 = rotor2.run_freq_response(speed_range)
# results1.plot_magnitude(inp = 4*4, out = 4*4)
# %%
