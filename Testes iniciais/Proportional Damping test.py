#%%
import numpy as np
import plotly.graph_objects as go
import ross as rs
steel = rs.materials.steel

L_total = 1
D = 0.002
n_el = 10

shaft = [
    rs.ShaftElement(
        1.0 / (n_el),
        idl=0,
        odl=D,
        material=steel,
        alpha = 8,
        beta = 1e-5,
        shear_effects=False,
        rotary_inertia=False,
    )
    for i in range(n_el)
]

bearings = [
    rs.BearingElement(0, kxx=1e15, cxx=0),
    rs.BearingElement(n_el, kxx=1e15, cxx=0),
]

rotor = rs.Rotor(shaft_elements=shaft, bearing_elements=bearings)
# %%
