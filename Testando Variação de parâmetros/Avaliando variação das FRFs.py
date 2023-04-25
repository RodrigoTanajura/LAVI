# Avaliando variação das FRFs conforme os parâmetros
#%%
import ross as rs
import numpy as np
from math import trunc
import time
from matplotlib import pyplot as plt

steel = rs.materials.steel

kxx = 8e5
kyy = 1e6

cxx = 12
cyy = 10

k = [kxx,kyy,cxx,cyy]

All_stats = []

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

speed_range = np.linspace(0,500,1001)

bounds=[(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

n = 5

# kxx
print("Varying kxx")
r = []
tempo = []
intervalo = np.linspace(bounds[0][0],bounds[0][1],n)

for i in range(len(intervalo)):
    st = time.time()
    bearings = [
        rs.BearingElement(0, kxx=intervalo[i], kyy=kyy,
        cxx=cxx, cyy=cyy),
        rs.BearingElement(6, kxx=intervalo[i], kyy=kyy,
        cxx=cxx, cyy=cyy)
    ]
    rotor1 = rs.Rotor(shaft, disks, bearings)
    results = rotor1.run_freq_response(speed_range)
    res= results.freq_resp[0][1*4]
    r.append(np.sqrt(res.imag**2 + res.real**2))
    et = time.time()
    tempo.append(et-st)
    print("Iteration {}, duration: {}".format(i+1, et-st))

for i in range(len(intervalo)):
    plt.plot(speed_range,r[i],label='={}'.format(round(intervalo[i],2)))
# plt.title('Average FRF calculation time: {}'.format(round(np.mean(tempo),2)))
plt.legend()
plt.show()
All_stats.append(r)
#%%
# kyy
print("Varying kyy")
r = []
tempo = []
intervalo = np.linspace(bounds[1][0],bounds[1][1],n)

for i in range(len(intervalo)):
    st = time.time()
    bearings = [
        rs.BearingElement(0, kxx=kxx, kyy=intervalo[i],
        cxx=cxx, cyy=cyy),
        rs.BearingElement(6, kxx=kxx, kyy=intervalo[i],
        cxx=cxx, cyy=cyy)
    ]
    rotor1 = rs.Rotor(shaft, disks, bearings)
    results = rotor1.run_freq_response(speed_range)
    res= results.freq_resp[0][1*4]
    r.append(np.sqrt(res.imag**2 + res.real**2))
    et = time.time()
    tempo.append(et-st)
    print("Iteration {}, duration: {}".format(i+1, et-st))

for i in range(len(intervalo)):
    plt.plot(speed_range,r[i],label='={}'.format(round(intervalo[i],2)))
plt.legend()
plt.show()
All_stats.append(r)
#%%
# cxx
print("Varying cxx")
r = []
tempo = []
intervalo = np.linspace(bounds[2][0],bounds[2][1],n)

for i in range(len(intervalo)):
    st = time.time()
    bearings = [
        rs.BearingElement(0, kxx=kxx, kyy=kyy,
        cxx=intervalo[i], cyy=cyy),
        rs.BearingElement(6, kxx=kxx, kyy=kyy,
        cxx=intervalo[i], cyy=cyy)
    ]
    rotor1 = rs.Rotor(shaft, disks, bearings)
    results = rotor1.run_freq_response(speed_range)
    res= results.freq_resp[0][1*4]
    r.append(np.sqrt(res.imag**2 + res.real**2))
    et = time.time()
    tempo.append(et-st)
    print("Iteration {}, duration: {}".format(i+1, et-st))

for i in range(len(intervalo)):
    plt.plot(speed_range,r[i],label='={}'.format(round(intervalo[i],2)))
plt.legend()
plt.show()
All_stats.append(r)
#%%
# cyy
print("Varying cyy")
r = []
tempo = []
intervalo = np.linspace(bounds[3][0],bounds[3][1],n)

for i in range(len(intervalo)):
    st = time.time()
    bearings = [
        rs.BearingElement(0, kxx=kxx, kyy=kyy,
        cxx=cxx, cyy=intervalo[i]),
        rs.BearingElement(6, kxx=kxx, kyy=kyy,
        cxx=cxx, cyy=intervalo[i])
    ]
    rotor1 = rs.Rotor(shaft, disks, bearings)
    results = rotor1.run_freq_response(speed_range)
    res= results.freq_resp[0][1*4]
    r.append(np.sqrt(res.imag**2 + res.real**2))
    et = time.time()
    tempo.append(et-st)
    print("Iteration {}, duration: {}".format(i+1, et-st))
for i in range(len(intervalo)):
    plt.plot(speed_range,r[i],label='={}'.format(round(intervalo[i],2)))
plt.legend()
plt.show()
All_stats.append(r)
# %%
