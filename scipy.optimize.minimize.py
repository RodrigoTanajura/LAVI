#%%
#Exemplo de aplicação de function minimization usando a biblioteca scipy
import numpy as np
from scipy.optimize import minimize

def objective_fcn(x):
    x1 = x[0]
    x2 = x[1]
    return x1**2 + x1*x2

def equality_constraint(x):
    x1 = x[0]
    x2 = x[1]
    return x1**3 + x1*x2 - 100
    
def inequality_constraint(x):
    x1 = x[0]
    x2 = x[1]
    return x1**2 +x2 - 50

bounds_x1 = (-100,100)
bounds_x2 = (-100,100)

bounds = [bounds_x1, bounds_x2]

constraint1 = {'type': 'eq', 'fun': equality_constraint}
constraint2 = {'type': 'ineq', 'fun': inequality_constraint}

constraint = [constraint1, constraint2]

x0 = [1,1]

result = minimize(objective_fcn,x0, method = 'SLSQP', bounds = bounds, constraints=constraint)

print(result)

# %%
