#%%
import ross as rs
import numpy as np
import plotly
from ross import FrequencyResponseResults
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from random import randint, gauss

shaft_elements = []
bearing_seal_elements = []
disk_elements = []
steel = rs.materials.steel

for i in range(6):
    shaft_elements.append(rs.ShaftElement(L=0.25, material=steel, n=i, idl=0, odl=0.05))

disk_elements.append(
    rs.DiskElement.from_geometry(n=2, material=steel, width=0.07, i_d=0.05, o_d=0.28)
)

disk_elements.append(
    rs.DiskElement.from_geometry(n=4, material=steel, width=0.07, i_d=0.05, o_d=0.35)
)
# bearing_seal_elements.append(rs.BearingElement(n=0, kxx=1e6, kyy=1e6, cxx=0, cyy=0))
# bearing_seal_elements.append(rs.BearingElement(n=6, kxx=1e6, kyy=1e6, cxx=0, cyy=0))

rotor591c = rs.Rotor(
    shaft_elements=shaft_elements,
    bearing_elements=bearing_seal_elements,
    disk_elements=disk_elements,
)

for i in range(1,6):
    print (i)
    bearing_seal_elements = []
    bearing_seal_elements.append(rs.BearingElement(n=0, kxx=i*1e6, kyy=i*1e6, cxx=0, cyy=0))
    bearing_seal_elements.append(rs.BearingElement(n=6, kxx=i*1e6, kyy=i*1e6, cxx=0, cyy=0))
    rotor591c = rs.Rotor(
    shaft_elements=shaft_elements,
    bearing_elements=bearing_seal_elements,
    disk_elements=disk_elements,
)
    rotor591c.plot_rotor()



# samples = 100
# speed_range = np.linspace(0, 500, samples) # rads/s
# results1 = rotor591c.run_freq_response(speed_range=speed_range)

# noise = []
# for i in range(samples):
#     noise.append(gauss(0,0.001))

# plt.scatter(speed_range, np.abs(results1.freq_resp[16][16]))

# def xp():
    

# parametros, pcov= curve_fit(, speed_range, np.abs(results1.freq_resp[16][16]), bounds = [0,np.inf])
# %%
