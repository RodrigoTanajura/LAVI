#%%
import ross as rs
import numpy as np

k = [855130.682384417,52024576.8035698,7.45161127921376,25.5867100796101,1197651.61618692,702261906.673702, 33.6793402012332,91.0335488825413,2.73016153382567, 4.84982845688298e-06,770.441680165024]
material = rs.Material("Aldemir", 7850, E=2.05e11, G_s=None, Poisson=0.29)
steel = rs.materials.steel

L=[0, 25, 49, 65, 80, 110, 140, 161, 190, 220, 250, 285, 295, 305, 330, 360, 390, 
420, 450, 480, 510, 523, 533, 543, 570, 594, 630, 664, 700, 730, 760, 790, 830, 862]
L = [L[i] - L[i - 1] for i in range(1, len(L))]

odl = 0.017 
list_odl = [odl]*len(L)
for i in [11,12,21,22]:
    list_odl[i] = list_odl[i]+0.020
list_odl[2] = 0.029
list_odl[3] = 0.029

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

disks = [rs.DiskElement(Id=0.007513248437500, Ip=0.003844540885417,
     m=2.6375, color='Firebrick', n=12),
    rs.DiskElement(Id=0.007547431937500, Ip=0.003862032635417,
     m=2.6495, color='Firebrick', n=22)
]

bearing0 = rs.BearingElement(
   n=3, kxx=k[0], kyy=k[1], cxx=k[2], cyy=k[3])
bearing1 = rs.BearingElement(
    n=30, kxx=k[4], kyy=k[5], cxx=k[6], cyy=k[7])


def width(massa):
    return (4*massa/(material.rho*np.pi*(odl/10)**2))

disks.append(rs.DiskElement.from_geometry(n=0, material=material, width=width(0.14945*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=1, material=material, width=width(0.14945*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=3, material=material, width=width(1.1252*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=30, material=material, width=width(0.0441*2), i_d=0, o_d=odl/10))

aldemir = rs.Rotor(shaft, disks, [bearing0,bearing1])
aldemir.save("aldemir.toml")
aldemir.plot_rotor()
# %%
