#%%
#Exemplo de aplicação de function minimization usando a biblioteca scipy
import numpy as np
import scipy.optimize as opt

x0 = [1,100]

def banana(x0):
    x = x0[0]
    y = x0[1]
    f = (1-x)**2 +100*(y-x**2)**2
    return f

optminimize = opt.minimize(banana, x0)
optfmin = opt.fmin(banana, x0)
#optminimizescalar = opt.minimize_scalar(banana(x0))
print(optminimize)
# %%