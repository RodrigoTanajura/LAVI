# Testando dados experimentais de uberlandia com DE
#%%
#%%

# Testando otimizacao em duas etapas. Primeiro buscamos encontrar as frequencias naturais do sistema (rigidezes) e 
# depois refazemos a otimizacao para ajustar a amplitude dos picos (amortecimento).
print("Iniciando bibliotecas e funções")
import time
st1 = time.time()
from matplotlib import pyplot as plt
import numpy as np
from numpy.linalg import norm 
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from math import trunc
from scipy.signal import find_peaks
import pandas as pd
import openpyxl
import scipy as sp
from scipy.optimize import minimize, least_squares
from scipy.optimize import differential_evolution, dual_annealing, shgo

arq1 = r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D2S2S4\TRAC7.TXT"

arquivos = [r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D1S1S3\TRAC9.TXT",
            r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D1S1S3\TRAC11.TXT",
            r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D1S2S4\TRAC1.TXT",
            r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D1S2S4\TRAC3.TXT",
            r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D2S1S3\TRAC13.TXT",
            r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D2S1S3\TRAC15.TXT",
            r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D2S2S4\TRAC5.TXT",
            r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Tentativa 2\TXT\D2S2S4\TRAC7.TXT"]

disco_1=12
disco_2=23
sensor_1=7
sensor_2=7
sensor_3=27
sensor_4=27

ordem = [(disco_1, sensor_1), (disco_1, sensor_3),(disco_1, sensor_2),(disco_1, sensor_4),
         (disco_2, sensor_1),(disco_2, sensor_3),(disco_2, sensor_2),(disco_2, sensor_4)]

pos = np.array(ordem)*4

for i in [0,1,4,5]:
    pos[i][1] += 1


def extract_data(arquivo):
    freq = np.loadtxt(arquivo[0])[:, 0]*2*np.pi
    all_z = []
    for i in arquivo:
        amp = np.loadtxt(i)[:, 1]
        im = np.loadtxt(i)[:, 2]
        z = np.sqrt(amp**2 + im**2)
        all_z.append(z[10:])
    return all_z, freq[10:]

amp, freq = extract_data(arquivos)

m1=0.79945 # massa do acoplamento na direcao x
m2=0.49045 # massa do acoplamento na direcao y
k1=1 # rigidez do acoplamento na direcao x
k2=1 # rigidez do acoplamento da direcao y

kxx = 1e6
kyy = 8e5
cxx = 5
cyy = 10

k = [kxx, kyy, k1, k2]

#Definindo o rotor de teste. Rotor simples de 6 nós com 1 disco central.
    
#                        A     A                 M1                      S1                      D1    
#                       0     1     2     3     4     5     6     7     8     9    10     11    12    13    14    15    16   17    18    19    20    21  
L = np.array([    0,   30,   60,   100,  120,  150,  180,  210,  244,  270,  300,  335,  345,  355,  375,  400,  430,  460,  490,  520,  550, 578,                     
                    616,  626,  636,   670,  700,  730,  765,  790,  820,  850,  880,  915,  950,  985, 1000])
#                       22    23    24     25    26    27    28    29    30    31    32    33   34     35    36      
#                             D2                             S2                            M2    
odl = 0.017
L = L*0.001
L = [L[i] - L[i - 1] for i in range(1, len(L))]

list_odl = [odl]*len(L)
for i in [11,12,22,23]:
    list_odl[i] += 0.011

Nnos = len(L)
Nele = Nnos-1

disco_1=12
disco_2=23
sensor_1=7
sensor_2=27

material = rs.Material("Aldemir", 7850, E=205e9,
                        G_s=None, Poisson=0.29, color='#525252')

shaft = []
for i in range(len(L)):
    shaft.append(
        rs.ShaftElement(
            material=material,
            L=L[i],
            idl=0,
            odl=list_odl[i],
            rotary_inertia=False,
            shear_effects=False
    ))

disks = [rs.DiskElement(Id=0.007415600636357646, Ip=0.003707800318178823,
    m=2.67429, color='Firebrick', n=disco_1),
    rs.DiskElement(Id=0.007377186232163936, Ip=0.003688593116081968,
    m=2.66029, color='Firebrick', n=disco_2)]

bounds=[(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]

bounds1 = [(trunc(0.8*1.1*kxx),trunc(1.2*1.1*kxx)),
        (trunc(0.8*1.1*kyy),trunc(1.2*1.1*kyy)),
        (0,150), (0,150)]

bounds2 = [(trunc(0.1*1.1*cxx),trunc(2*1.1*cxx)),
        (trunc(0.1*1.1*cyy),trunc(2*1.1*cyy))]
# mudar bounds2
# X0 é o vetor que marca o ponto médio das condicoes de contorno.

x0 = ((bounds[0][0]+bounds[0][1])/2,(bounds[1][0]+bounds[1][1])/2,
(bounds[2][0]+bounds[2][1])/2, (bounds[3][0]+bounds[3][1])/2 )

def FRF_only_stiffness(k, speed_range, pos):
    FRF_resp = np.zeros(len(speed_range))
    bearings = [
            rs.BearingElement(4, kxx=k[0], kyy=k[1],
            cxx=x0[2], cyy=x0[3]),
            rs.BearingElement(33, kxx=k[0], kyy=k[1],
            cxx=x0[2], cyy=x0[3])
        ]
    rotor = rs.Rotor(shaft, disks, bearings)

    M_mod = rotor.M() # Alterando a matriz para contabilizar o acoplamento
    M_mod[0*4,0*4] += m1 #horizontal
    M_mod[0*4+1,0*4+1] += m2 #vertical
    M_mod[1*4,1*4] += m1 #horizontal
    M_mod[1*4+1,1*4+1] += m2 #vertical
    
    K_mod = np.zeros(M_mod.shape)
    K_mod[0*4+2,0*4+2] += k[2]
    K_mod[0*4+3,0*4+3] += k[3]
    K_mod[1*4+2,1*4+2] += k[2]
    K_mod[1*4+3,1*4+3] += k[3]
    K0 = rotor.K(0) + K_mod
    Lambda, Phi = sp.linalg.eigh(K0, M_mod)
    Phi = Phi[:,:10]
    M_red = Phi.T@M_mod@Phi
    G_red = Phi.T@rotor.G()@Phi
    for i in range(len(speed_range)):
        w=speed_range[i]
        K_red = Phi.T @ (rotor.K(w) + K_mod) @ Phi
        C_red = Phi.T @ rotor.C(w) @ Phi
        H_red = -w**2*M_red + (1j*w)*(C_red+w*G_red) + K_red
        Hinv_red = sp.linalg.inv(H_red)
        Hinv = Phi@Hinv_red@Phi.T
        FRF_resp[i] = np.abs(Hinv[pos[0]][pos[1]])
    return FRF_resp

def FRF_only_dampening(c, speed_range, pos):
    FRF_resp = np.zeros(len(speed_range))
    bearings = [
            rs.BearingElement(4, kxx=x1[0], kyy=x1[1],
            cxx=c[0], cyy=c[1]),
            rs.BearingElement(33, kxx=x1[0], kyy=x1[1],
            cxx=c[0], cyy=c[1])
        ]
    rotor = rs.Rotor(shaft, disks, bearings)

    M_mod = rotor.M() # Alterando a matriz para contabilizar o acoplamento
    M_mod[0*4,0*4] += m1 #horizontal
    M_mod[0*4+1,0*4+1] += m2 #vertical
    M_mod[1*4,1*4] += m1 #horizontal
    M_mod[1*4+1,1*4+1] += m2 #vertical
    K_mod = np.zeros(M_mod.shape)
    K_mod[0*4+2,0*4+2] += x1[2]
    K_mod[0*4+3,0*4+3] += x1[3]
    K_mod[1*4+2,1*4+2] += x1[2]
    K_mod[1*4+3,1*4+3] += x1[3]
    K0 = rotor.K(0) + K_mod

    Lambda, Phi = sp.linalg.eigh(K0, M_mod)
    Phi = Phi[:,:10]
    M_red = Phi.T@M_mod@Phi
    G_red = Phi.T@rotor.G()@Phi
    for i in range(len(speed_range)):
        w=speed_range[i]
        K_red = Phi.T @ rotor.K(w) @ Phi
        C_red = Phi.T @ rotor.C(w) @ Phi

        H_red = -w**2*M_red + (1j*w)*(C_red+w*G_red) + K_red
        Hinv_red = sp.linalg.inv(H_red)
        Hinv = Phi@Hinv_red@Phi.T
        FRF_resp[i] = np.abs(Hinv[pos[0]][pos[1]])
    return FRF_resp

# Localizando picos e alterando o speed range pra focar neles.

def peak_finder(amp, freq, intervalo=10):
    speed_range = []
    amp2 = []
    for j in amp:
        picos = find_peaks(j, prominence=5e-6)
        speed = freq
        for i in range(len(picos[0])):
            speed_range = np.append(speed_range, speed[picos[0][i]-intervalo:picos[0][i]+intervalo])
            amp2 = np.append(amp2, amp[picos[0][i]-intervalo:picos[0][i]+intervalo])
    return amp2, speed_range

z1, speed_range = peak_finder(amp, freq)

# Funcao objetivo Ritto
def objective_stiffness(k, speed_range, pos):
    res = []
    for i in range(len(pos)):
        res.append((norm(z1[i] - FRF_only_stiffness(k, speed_range[i], pos[i], 2)**2)/(norm(z1[i], 2)**2)))
    return sum(res)

def objective_dampening(c, speed_range, pos):
    res = (norm(z1 - FRF_only_dampening(c, speed_range, pos), 2)**2)/(norm(z1, 2)**2)
    return res
# return da funcao objetivo para amortecimento é multiplicado por 100, para aumentar a sensibilidade da resposta. Ou seja, a funcao objetivo terá de ser minimizada 100 vezes mais que para a rigidez.

# # Otimizações:
# # Somente as rigidezes:
# print("Otimizando a rigidez")
# st1 = time.time()
# res = differential_evolution(objective_stiffness, bounds1, disp=True, polish=False, args=(speed_range, pos))
# et1 = time.time()

# x1 = [res.x[0], res.x[1], res.x[2], res.x[3]] # Fixando as rigidezes com o resultado da otimizacao acima. Esses dados serao usados a seguir na otimizacao do amortecimento. 

# print("Otimizando amortecimento")
# # Agora somente o amortecimento:

# res_dampening = differential_evolution(objective_dampening, bounds2, disp=True, polish=True, args=(speed_range, pos))
# et2 = time.time()

# print('Tempo total:{}\n Tempo 1:{}\n Tempo 2:{}\n'.format(-st1+et2,et1-st1,et2-et1) )

# print("Valores de rigidez:{}\n".format(res.x))
# print("Valores de amortecimento:{}".format(res_dampening.x))

# %%
