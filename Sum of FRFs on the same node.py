# Sum of FRFs on the same node
#%%
import numpy as np

frfs = [np.arange(0,100),
        np.arange(100,200),
        np.arange(200,300),
        np.arange(300,400)]

def sum_FRF(list_of_frf):
    return sum(frfs[i] for i in range(len(frfs)))/len(frfs)
# %%
