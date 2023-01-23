#%%
#Obs: aplicar clustering points para obter maior densidade de pontos em áreas mais 
# relevantes.
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from scipy.optimize import curve_fit, minimize_scalar, minimize
   

steel = rs.materials.steel
samples = 1001
speed_range = np.linspace(0, 150, samples) # rads/s
k = 1e6

def create_rotor(kxx):
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
        rs.BearingElement(0, kxx=kxx, cxx=0),
        rs.BearingElement(6, kxx=kxx, cxx=0),
    ]
    rotor = rs.Rotor(shaft, disks, bearings)
    results1 = rotor.run_freq_response(speed_range)
    return np.abs(results1.freq_resp[16][16])

noise = []
for i in range(samples):
    noise.append(gauss(0,0.0000))

def real_input(kxx):
    return create_rotor(kxx) + noise

def objective(kxx):
    res = create_rotor(kxx) - real_input(k)
    return np.linalg.norm(res)

result = minimize_scalar(objective)
p_at_minimum = result.x  # i'm a scalar!
print('minimize_scalar result: ', p_at_minimum)



plt.plot(speed_range, create_rotor(k),'b')
plt.plot(speed_range, create_rotor(78),'r')
plt.plot(speed_range, real_input(k),'o')
plt.show()
# x0 = [1]
# res = minimize(opt, x0, method='nelder-mead',
#                options={'xatol': 1e-8, 'disp': True})

# parametros, pcov= curve_fit(create_rotor(speed_range,k), speed_range,  opt(speed_range,k), bounds = [0,np.inf])


# plt.scatter(speed_range, np.abs(create_rotor(1e6))+noise)
# %%
