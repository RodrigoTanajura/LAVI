#%%
from matplotlib import pyplot as plt
import numpy as np
import plotly.graph_objects as go
import ross as rs
import plotly.io as pio
from random import randint, gauss
from math import trunc
from PIL import Image

arq1 = r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\TXT\0 Hz\Tentativa 2\D2S2S4\TRAC7.TXT"
arq2 = r"C:\Users\diwiw\OneDrive\Documentos\A UFRJ\Estágio - LAVI\Dados de Uberlandia\Raimundo\FRF_VD2\TRAC2_D2_P4_V.TXT"
title = "FRF w = 0."

def plotfile1(arquivo):
    freq = np.loadtxt(arquivo)[:, 0]
    amp = np.loadtxt(arquivo)[:, 1]
    im = np.loadtxt(arquivo)[:, 2]
    z = np.sqrt(amp**2 + im**2)
    return [freq[10:], z[10:]]

a1 = plotfile1(arq1)
def plotfile2(arquivo):
    freq = np.loadtxt(arquivo)[:, 0]
    amp = np.loadtxt(arquivo)[:, 1]
    im = np.loadtxt(arquivo)[:, 2]
    z = np.sqrt(amp**2 + im**2)
    return [freq[10:], z[10:]*1e-3]
a2 = plotfile2(arq2)

plt.semilogy(a1[0], a1[1])
# plt.semilogy(a2[0], a2[1], label="Raimundo")
plt.grid()
plt.xlabel("Hertz")
plt.ylabel("Metros")
plt.legend()
plt.title(title)
plt.savefig(title)
# %%
