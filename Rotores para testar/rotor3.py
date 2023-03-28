#%%
import ross as rs
import numpy as np
from scipy import signal

steel = rs.materials.steel

kxx = 1.2e6
kyy = 1e6

cxx = 80
cyy = 60


k = [kxx,kyy,cxx,cyy]

shaft = [
        rs.ShaftElement(1.2 / (12), idl=0.03, odl=0.08, material=steel)
        for i in range(12)
    ]

disks = [
    rs.DiskElement.from_geometry(
        n=(0), material=steel, width=0.08, i_d=0.08, o_d=0.4
    ),
    rs.DiskElement.from_geometry(
        n=(1), material=steel, width=0.04, i_d=0.08, o_d=0.3
    ),
    rs.DiskElement.from_geometry(
        n=(2), material=steel, width=0.02, i_d=0.08, o_d=0.2
    )
]

bearings = [
        rs.BearingElement(12, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3]),
        rs.BearingElement(3, kxx=k[0], kyy=k[1],
        cxx=k[2], cyy=k[3])
    ]

speed_range = np.linspace(0,500, 501)

rotor3 = rs.Rotor(shaft, disks, bearings)
rotor3.plot_rotor()
results1 = rotor3.run_freq_response(speed_range)
results1.plot_magnitude(inp = 0, out = 1*4)
# %%
