#%%
import ross as rs

steel = rs.materials.steel

kxx = 8e5
kyy = 1e6

cxx = 12
cyy = 10

k = [kxx,kyy,cxx,cyy]


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
