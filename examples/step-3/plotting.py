# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # H5 files and plotting

# %%
import numpy as np
import h5py

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from scipy.interpolate import griddata
from matplotlib import cm

# %%
refinement = 5
filename = "solution_%d.h5" % refinement
file = h5py.File(filename, "r")

# %%
for key, value in file.items():
    print(key, " : ", value)

# %%
nodes = np.array(file["/nodes"])
cells = np.array(file["/cells"])
solution = np.array(file["/solution"])

x, y = nodes.T

# %%
nodes

# %%
cells

# %%
solution

# %%
x

# %%
y

# %%
cell_x = x[cells.flatten()]
cell_y = y[cells.flatten()]

# %%
cell_x

# %%
cell_y

# %%
cell_ids = np.repeat(np.arange(cells.shape[0]), 4)
cell_ids

# %%
n_cells = cells.shape[0]
n_cells

# %%
meshdata = pd.DataFrame({"x": cell_x, "y": cell_y, "ids": cell_ids})

# %%
meshdata

# %%
fig, ax = plt.subplots()
ax.set_aspect("equal", "box")
ax.set_title("grid at refinement level #%d" % refinement)

for cell_id, cell in meshdata.groupby(["ids"]):
    cell = pd.concat([cell, cell.head(1)])
    plt.plot(cell["x"], cell["y"], c="k")

# %% [markdown]
# An alternative is to use the `matplotlib` built-in `Polygon` class

# %%
fig, ax = plt.subplots()
ax.set_aspect("equal", "box")
ax.set_title("grid at refinement level #%d" % refinement)
for cell_id, cell in meshdata.groupby(["ids"]):
    patch = Polygon(cell[["x", "y"]], facecolor="w", edgecolor="k")
    ax.add_patch(patch)


# %% [markdown]
# ## A color plot of the solution

# %%
nx = int(np.sqrt(n_cells)) + 1
nx *= 10
xg = np.linspace(x.min(), x.max(), nx)
yg = np.linspace(y.min(), y.max(), nx)

xgrid, ygrid = np.meshgrid(xg, yg)
solution_grid = griddata((x, y), solution.flatten(), (xgrid, ygrid), method="linear")

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.set_title("solution at refinement level #%d" % refinement)
c = ax.pcolor(xgrid, ygrid, solution_grid, cmap=cm.viridis)
fig.colorbar(c, ax=ax)

plt.show()

# %% [markdown]
# ## Convergence

# %%
mean_values = np.zeros((8,))
point_values = np.zeros((8,))
dofs = np.zeros((8,))

for refinement in range(1, 9):
    filename = "solution_%d.h5" % refinement
    file = h5py.File(filename, "r")
    mean_values[refinement - 1] = np.array(file["/mean_value"])[0]
    point_values[refinement - 1] = np.array(file["/point_value"])[0]
    dofs[refinement - 1] = np.array(file["/solution"]).shape[0]


# %%
mean_values

# %%
point_values

# %%
dofs

# %%
mean_error = np.abs(mean_values[1:] - mean_values[:-1])
plt.loglog(dofs[:-1], mean_error)
plt.grid()
plt.xlabel("#DoFs")
plt.ylabel("mean value error")
plt.show()

# %%
point_error = np.abs(point_values[1:] - point_values[:-1])
plt.loglog(dofs[:-1], point_error)
plt.grid()
plt.xlabel("#DoFs")
plt.ylabel("point value error")
plt.show()

# %% [markdown]
# # Using *python* equivalent of *ggplot* of R
#
# **[plotnine](https://plotnine.readthedocs.io/en/v0.12.4/)**

# %%
# !pip install plotnine

# %%
from plotnine import (
    ggplot,
    aes,
    geom_raster,
    geom_polygon,
    geom_line,
    labs,
    scale_x_log10,
    scale_y_log10,
    ggtitle,
)  # noqa: E402

# %%
plot = (
    ggplot(meshdata, aes(x="x", y="y", group="ids"))
    + geom_polygon(fill="white", colour="black")
    + ggtitle("grid at refinement level #%d" % refinement)
)

print(plot)

# %%
colordata = pd.DataFrame({"x": x, "y": y, "solution": solution.flatten()})
colordata

# %%
plot = (
    ggplot(colordata, aes(x="x", y="y", fill="solution"))
    + geom_raster(interpolate=True)
    + ggtitle("solution at refinement level #%d" % refinement)
)

print(plot)

# %%
convdata = pd.DataFrame(
    {"dofs": dofs[:-1], "mean_value": mean_error, "point_value": point_error}
)

convdata

# %%
plot = (
    ggplot(convdata, mapping=aes(x="dofs", y="mean_value"))
    + geom_line()
    + labs(x="#DoFs", y="mean value error")
    + scale_x_log10()
    + scale_y_log10()
)

plot.save("mean_error.pdf", dpi=60)
print(plot)

# %%
plot = (
    ggplot(convdata, mapping=aes(x="dofs", y="point_value"))
    + geom_line()
    + labs(x="#DoFs", y="point value error")
    + scale_x_log10()
    + scale_y_log10()
)

plot.save("point_error.pdf", dpi=60)
print(plot)
