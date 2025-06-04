import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys
import plotstyle

import BasisFunctions

from matplotlib.colors import hsv_to_rgb

# Ensure all input files are correctly referenced
try:
    run_number = sys.argv[1]
    beta = sys.argv[2]
except:
    run_number = input('Run number: ')
    beta = input('Beta: ')

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# Function to get data
def get_data(filename, variables):
    df = pd.read_csv(filename, delim_whitespace=True, engine='python', names=variables)
    return df

# Formatting run number
run_number = f"Run{int(run_number):03d}"
Rvecsfile = os.path.join(run_number, "Rvecs.txt")
configfile = os.path.join(run_number, f"configs_{beta}.txt")

# Load Rvecs and config data
data = get_data(Rvecsfile, ["x", "y", "z"])
Rvecs = np.array(data, dtype=float)

config = []
with open(configfile, 'r') as infile:
    for line in infile:
        config.append(np.array(line.split(), dtype=float))

config = np.concatenate(config).reshape(-1, 3)

#config = np.zeros((Rvecs.shape[0], 3))
#config[:, 2] = -np.ones(Rvecs.shape[0])

# Normalize and calculate colors outside the loop
norm = plt.Normalize(vmin=-1, vmax=1)
colors = plt.cm.viridis(norm(config[:, 2]))

fig, axs = plt.subplots(1, 3, figsize=(12, 5))

for ax in axs:
    ax.set_aspect('equal')

# Calculate marker size
x_min, x_max = np.min(Rvecs[:, 0]), np.max(Rvecs[:, 0])
y_min, y_max = np.min(Rvecs[:, 1]), np.max(Rvecs[:, 1])
x_range = x_max - x_min
y_range = y_max - y_min
#marker_size = min(fig.get_size_inches()[0] * fig.dpi / x_range, fig.get_size_inches()[1] * fig.dpi / y_range) * 1.5/100


# Generate the grid and color matrix
x_unique = np.unique(Rvecs[:, 0])
y_unique = np.unique(Rvecs[:, 1])
x_indices = np.searchsorted(x_unique, Rvecs[:, 0])
y_indices = np.searchsorted(y_unique, Rvecs[:, 1])

Nx = len(x_unique)
Ny = len(y_unique)

# Find nbrs in x and y direction:
nbrs = [[BasisFunctions.coord_to_site(Rvecs[i]+np.array([1,0,0]), Nx, Ny), BasisFunctions.coord_to_site(Rvecs[i]+np.array([0,1,0]), Nx, Ny)] for i in range(len(Rvecs))]

def colors(conf1, conf2, conf3):
    # Calculate hue from xy components
    xy_angles = np.arctan2(conf2, conf1)
    hues = (xy_angles + np.pi) / (2 * np.pi)

    # Calculate lightness from z component
    lightness = (conf3 + 1) / 2  # Normalize z to [0, 1]; z = -1 maps to 0 (black), z = 1 maps to 1 (white)

    # Create colors using HSV, with full saturation
    colors = np.zeros((config.shape[0], 3))  # RGB
    for i in range(len(hues)):
        h = hues[i]
        z = conf3[i]
        if z >= 0:
            # Interpolate linearly between the HSV color and white
            rgb = hsv_to_rgb([h, 1, 1])
            colors[i, :3] = rgb * (1 - z) + np.array([1, 1, 1]) * z
        else:
            # Interpolate linearly between the HSV color and black
            rgb = hsv_to_rgb([h, 1, 1])
            colors[i, :3] = rgb * (1 + z)
        #colors[i, 3] = 1  # Fully opaque

    color_grid = np.zeros((len(y_unique), len(x_unique), 3))

    for color, x_idx, y_idx in zip(colors, x_indices, y_indices):
        color_grid[y_idx, x_idx] = color

    # Reverse the y-axis to match the typical matrix representation
    color_grid = np.flip(color_grid, axis=0)

    return color_grid

marker_size = 5

xycolor_grid = colors(config[:, 0], config[:, 1], config[:, 2])
yzcolor_grid = colors(config[:, 1], config[:, 2], config[:, 0])
zxcolor_grid = colors(config[:, 2], config[:, 0], config[:, 1])

X, Y = np.meshgrid(x_unique, y_unique)

# Use scatter function directly with arrays
axs[0].imshow(xycolor_grid, aspect='equal', extent=(x_unique.min() - 0.5, x_unique.max() + 0.5, y_unique.min() - 0.5, y_unique.max() + 0.5))
axs[1].imshow(yzcolor_grid, aspect='equal', extent=(x_unique.min() - 0.5, x_unique.max() + 0.5, y_unique.min() - 0.5, y_unique.max() + 0.5))
axs[2].imshow(zxcolor_grid, aspect='equal', extent=(x_unique.min() - 0.5, x_unique.max() + 0.5, y_unique.min() - 0.5, y_unique.max() + 0.5))

for ax, t in zip(axs, ["xy", "yz", "zx"]):
    ax.set_xlim(x_min - 0.5, x_max + 0.5)
    ax.set_ylim(y_min - 0.5, y_max + 0.5)
    ax.grid(False)
    ax.set_title(f"{t}-plane")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

print(np.mean(abs(config[:,2])))

plt.show()