# -*- coding: utf-8 -*-
"""
Script to plot the results from example2_pbe_2d_fv.f90

"""

# %% Dependencies

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.animation import PillowWriter
import pandas as pd

# %% Import numerical results

relpath_output = '../output/example2'
filenames = {
    "u": 'u.txt',
    "x1": 'x1.txt',
    "x2": 'x2.txt'
}

data = {}
here = os.path.dirname(os.path.abspath(__file__))
abspath_output = os.path.join(here, relpath_output)
for var, filename in filenames.items():
    data[var] = pd.read_csv(os.path.join(abspath_output, filename),
                            delim_whitespace=True)

# Extract arrays
x1 = data['x1']['x(i)'].values
x2 = data['x2']['x(i)'].values
t = data['u']['t'].values
u = data['u'].iloc[:, 1:].values

tpoints = u.shape[0]
upoints = u.shape[1]

# %% Make animated surface plot

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
X1, X2 = np.meshgrid(x1, x2, indexing='ij')

text = ax.text(0.8*np.max(x1), np.max(x2), 1.5*np.max(u), "",
               bbox=dict(facecolor='white', alpha=0.8))

ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('u(x,t)')
ax.set_title("Solution of 2D PBE with 5th order WENO")

ax.set_xlim(round(np.min(x1)), round(np.max(x1)))
ax.set_ylim(round(np.min(x2)), round(np.max(x2)))
ytol = 0.1
ax.set_zlim(0, np.max(u)+ytol)
ax.grid(True)

fps = 20
writer = PillowWriter(fps=fps)
with writer.saving(fig, os.path.join(abspath_output, "example2d.gif"), 150):
    for i, ti in enumerate(t):
        ax.plot_surface(X1, X2, np.reshape(u[i, :], X1.shape, order='F'),
                        cmap=cm.coolwarm,
                        linewidth=0, antialiased=False)
        text.set_text(f"time = {ti:.2f}")
        writer.grab_frame()
        if i == (tpoints-1):
            for ii in range(2*fps):
                writer.grab_frame()
