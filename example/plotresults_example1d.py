# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# %% Dependencies

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from matplotlib.animation import PillowWriter
import pandas as pd

# %% Import results

filenames = {
    "u": "../output/u.txt",
    "xgrid": "../output/xgrid.txt"
}

data = {}
for var, filename in filenames.items():
    data[var] = pd.read_csv(filename, delim_whitespace=True)

# Extract arrays
x = data['xgrid']['x(i)'].values
t = data['u']['t'].values
u = data['u'].iloc[:, 1:].values

# %% Make surface plot

# Plot the surface
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
T, X = np.meshgrid(t, x)
surf = ax.plot_surface(X, T, u.T, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
# ax.view_init(35, 215)

# Customize the z axis.
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u(x,t)')
ax.zaxis.set_major_locator(LinearLocator(5))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=10, location='left')

# %% Make moving 2d plot

fig, ax = plt.subplots(1, 1)
line, = ax.plot([], [])
text = ax.text(2.5, 0.75, "", bbox=dict(facecolor='white', alpha=0.8))
ax.set_xlabel('x')
ax.set_ylabel('u(x,t)')
ax.set_title("Solution of Burger's equation with 5th order WENO")
ax.set_xlim(round(np.min(x)), round(np.max(x)))
ytol = 1.05
ax.set_ylim(ytol*round(np.min(u[0, :])), ytol*round(np.max(u[0, :])))
ax.grid(True)

writer = PillowWriter(fps=15)
with writer.saving(fig, "../output/example1d.gif", 100):
    for i, ti in enumerate(t):
        line.set_data(x, u[i, :])
        text.set_text(f"time = {ti:.2f}")
        writer.grab_frame()
