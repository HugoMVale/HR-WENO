# -*- coding: utf-8 -*-
"""
Script to plot the results from example1_burgers_1d_fv.f90

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

relpath_output = '../output/example1'
filenames = {
    "u": 'u.txt',
    "xgrid": 'x.txt'
}

data = {}
here = os.path.dirname(os.path.abspath(__file__))
abspath_output = os.path.join(here, relpath_output)
for var, filename in filenames.items():
    data[var] = pd.read_csv(os.path.join(abspath_output, filename),
                            delim_whitespace=True)

# Extract arrays
x = data['xgrid']['x(i)'].values
t = data['u']['t'].values
u = data['u'].iloc[:, 1:].values

tpoints = u.shape[0]
xpoints = u.shape[1]

# %% Make surface plot

title = "Solution of Burger's equation with 5th order WENO"

# Plot the surface
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
T, X = np.meshgrid(t, x)
surf = ax.plot_surface(X, T, u.T, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u(x,t)')
ax.set_title(title)
# ax.zaxis.set_major_locator(LinearLocator(3))
fig.colorbar(surf, shrink=0.5, aspect=10, location='left')

# %% Make animated line plot

fig, ax = plt.subplots(1, 1)
line, = ax.plot([], [], '-o', fillstyle='none')
text = ax.text(2.5, 0.75, "", bbox=dict(facecolor='white', alpha=0.8))
ax.set_xlabel('x')
ax.set_ylabel('u(x,t)')
ax.set_title(title)
ax.set_xlim(round(np.min(x)), round(np.max(x)))
ytol = 0.1
ax.set_ylim(np.min(u[0, :])-ytol, np.max(u[0, :])+ytol)
ax.grid(True)

fps = 20
writer = PillowWriter(fps=fps)
with writer.saving(fig, os.path.join(abspath_output, "example1d.gif"), 150):
    for i, ti in enumerate(t):
        line.set_data(x, u[i, :])
        text.set_text(f"time = {ti:.2f}")
        writer.grab_frame()
        if i == (tpoints-1):
            for ii in range(2*fps):
                writer.grab_frame()
