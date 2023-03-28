#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from math import trunc

# O objetivo desse script é simular o rotor do Aldemir e obter resultados próximos
# aos da primeira FRF.

#Dados da primeira FRF. Impacto no nó 13, sensor no nó 8. Ambos na direcao x.
freq = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 0]
amp = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 1]
im = np.loadtxt(r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\TESE Aldemir\TRAC37.TXT")[:, 2]

odl = 0.017 
z = np.sqrt(amp**2 + im**2)
z = z[25:,]

speed_range = []
z1 = []

# Fator de correcao de Hz para radianos.
j = 0 
f_de_correcao = np.pi*2
freq = freq[25:,]* f_de_correcao

# A variavel reducao_de_densidade funciona para diminuir a quantidade de pontos que o ROSS ira rodar.

reducao_de_densidade = 8
for i in range(int(len(freq)/reducao_de_densidade)):
    speed_range.append(freq[j])
    z1.append(z[j])
    j += reducao_de_densidade

# speed_range = speed_range[0:100]
# z1 = z1[0:100]

#Dados otimizados do rotor. Seriam o "gabarito".
k = [855130.682384417,52024576.8035698,7.45161127921376,25.5867100796101,1197651.61618692,702261906.673702, 33.6793402012332,91.0335488825413,2.73016153382567, 4.84982845688298e-06,770.441680165024]

# Outros possíveis valores para testar:
# k = [5e3,5e3,30,30,5e5,5e5,50,50]
# k = [5e6,5e7,1e1,1e1,5e9,5e9,1e1,1e1]
# k = [10,10,2.5,2.5,10,10,2.5,2.5]

#Instanciando o rotor.
material = rs.Material("Aldemir", 7850, E=2.05e11, G_s=None, Poisson=0.29)
steel = rs.materials.steel

L=[0, 25, 49, 65, 80, 110, 140, 161, 190, 220, 250, 285, 295, 305, 330, 360, 390, 
420, 450, 480, 510, 523, 533, 543, 570, 594, 630, 664, 700, 730, 760, 790, 830, 862]
L = [L[i] - L[i - 1] for i in range(1, len(L))]

# Alteracoes de diametro ao longo do eixo:
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

# Valores de massa, Id e Ip tirados diretamente do arquivo FEmodel.m linha 124.
disks = [rs.DiskElement(Id=0.007513248437500, Ip=0.003844540885417,
     m=2.6375, color='Firebrick', n=12),
    rs.DiskElement(Id=0.007547431937500, Ip=0.003862032635417,
     m=2.6495, color='Firebrick', n=22)
]

bearing0 = rs.BearingElement(
   n=3, kxx=k[0], kyy=k[1], cxx=k[2], cyy=k[3])
bearing1 = rs.BearingElement(
    n=30, kxx=k[4], kyy=k[5], cxx=k[6], cyy=k[7])
bearing2 = rs.BearingElement(
    n=0, kxx=k[10], kyy=k[10], cxx=0, cyy=0)
# A funcao width garante que o disco tenha a massa correta. Ela calcula o width necessário para 
# que o disco tenha a massa correta.

def width(massa):
    return (4*massa/(material.rho*np.pi*(odl/10)**2))

# Usei discos para simular as 4 massas pontuais. Defini o diametro como 1 decimo do diametro do eixo.

disks.append(rs.DiskElement.from_geometry(n=0, material=material, width=width(0.14945*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=1, material=material, width=width(0.14945*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=3, material=material, width=width(1.1252*2), i_d=0, o_d=odl/10))
disks.append(rs.DiskElement.from_geometry(n=30, material=material, width=width(0.0441*2), i_d=0, o_d=odl/10))

aldemir = rs.Rotor(shaft, disks, [bearing0,bearing1,bearing2])
aldemir.save("aldemir.toml")
aldemir.plot_rotor()

#Plotando os resultados:
results1 = aldemir.run_freq_response(speed_range)
results1.plot_magnitude(inp = 13*4, out = 8*4)
# %%
