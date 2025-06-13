import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import seaborn as sns
import os
import sys
import BasisFunctions
#import plotstyle

from matplotlib.colors import hsv_to_rgb


# Function to get data
def get_data(filename, variables):
    df = pd.read_csv(filename, delim_whitespace=True, engine='python', names=variables)
    return df

def colors(conf1, conf2, conf3):
    # Calculate hue from xy components
    xy_angles = np.arctan2(conf2, conf1)
    hues = (xy_angles + np.pi) / (2 * np.pi)

    # Create colors using HSV, with full saturation
    colors = np.zeros((len(conf1), 3))  # RGB
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

    try:
        color_grid = np.zeros((len(y_unique), len(x_unique), 3))
        for color, x_idx, y_idx in zip(colors, x_indices, y_indices):
            color_grid[y_idx, x_idx] = color
    except:
        n = int(np.sqrt(len(conf1)))
        print(n)
        color_grid = np.zeros((n, n, 3))
        y = 0
        x = 0
        for color in colors:
            color_grid[y, x] = color
            x +=1
            if x == n:
                x = 0
                y +=1

    # Reverse the y-axis to match the typical matrix representation
    color_grid = np.flip(color_grid, axis=0)

    try: 
        #Transate the matrix
        color_grid = np.roll(color_grid, tx, axis=1)
        color_grid = np.roll(color_grid, -ty, axis=0)
    except:
        pass

    return color_grid

def RingCondition(Q, theta):
    J1 = -1
    J2 = 0.28
    J3 = J2/2.
    qx = Q*np.cos(theta); qy = Q*np.sin(theta);
    return np.cos(qx) + np.cos(qy) - (abs(J1)/(2.0*J2))

def Q(theta):
    Qa, Qb = 0, np.pi    
    while True:
        Qc = (Qa + Qb) / 2.0
        if abs(RingCondition(Qc, theta)) < 1e-3:
            break
        if RingCondition(Qa, theta) * RingCondition(Qc, theta) > 0:
            Qa = Qc
        else:
            Qb = Qc
    return np.array([Qc * np.cos(theta), Qc * np.sin(theta)])


if __name__ == "__main__":

    # Ensure all input files are correctly referenced
    try:
        run_number = sys.argv[1]
        betas = sys.argv[2]
        tx = int(sys.argv[-2])
        ty = int(sys.argv[-1])
    except:
        run_number = input('Run number: ')
        betas = input('Beta: ')

    import matplotlib
    matplotlib.rc('xtick', labelsize=14)
    matplotlib.rc('ytick', labelsize=14)
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'

    # Formatting run number
    run_number = f"Run{int(run_number):03d}"
    Rvecsfile = os.path.join(run_number, "Rvecs.txt")
    configfile = os.path.join(run_number, "equconfigs.txt")

    # Load Rvecs and config data
    data = get_data(Rvecsfile, ["x", "y", "z"])
    Rvecs = np.array(data, dtype=float)

    config = []
    with open(configfile, 'r') as infile:
        lines = infile.readlines()

    theta1 = np.linspace(0, 2*np.pi, 1000)
    Qs1 = np.zeros((len(theta1), 2))
    for i in range(len(theta1)):
        Qs1[i,:] = Q(theta1[i])

    for i in range(0, len(lines), 2):
    #for i in range(len(lines)-2, len(lines), 2):
            beta = lines[i].strip()  # Fjern eventuelle ledende/trailende mellomrom

            if True: #beta in betas:
                config = np.array( lines[i+1].split(), dtype=float)
                config = config.reshape(-1, 3)

                # Normalize and calculate colors outside the loop
                norm = plt.Normalize(vmin=-1, vmax=1)

                fig, axs = plt.subplots(2, 4, figsize=(10, 6))
                fig.subplots_adjust(wspace=0.4)

                for ax in np.ndarray.flatten(axs):
                    ax.set_aspect('equal')

                # Generate the grid and color matrix
                x_unique = np.unique(Rvecs[:, 0])
                y_unique = np.unique(Rvecs[:, 1])
                x_indices = np.searchsorted(x_unique, Rvecs[:, 0])
                y_indices = np.searchsorted(y_unique, Rvecs[:, 1])

                Nx = len(x_unique)
                Ny = len(y_unique)

                # Find nbrs in x and y direction:
                nbrs = [[BasisFunctions.coord_to_site(Rvecs[i]+np.array([1,0,0]), Nx, Ny), BasisFunctions.coord_to_site(Rvecs[i]+np.array([0,1,0]), Nx, Ny)] for i in range(len(Rvecs))]

                xycolor_grid = colors(config[:, 0], config[:, 1], config[:, 2])
                yzcolor_grid = colors(config[:, 1], config[:, 2], config[:, 0])
                zxcolor_grid = colors(config[:, 2], config[:, 0], config[:, 1])

                X, Y = np.meshgrid(x_unique, y_unique)

                x_min, x_max = x_unique.min(), x_unique.max()
                y_min, y_max = y_unique.min(), y_unique.max()

                # Use scatter function directly with arrays
                axs[0,0].imshow(xycolor_grid, aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))
                axs[0,1].imshow(yzcolor_grid, aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))
                axs[0,2].imshow(zxcolor_grid, aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))

                for ax, t in zip(axs[0,:3], ["xy", "yz", "zx"]):
                    ax.set_xlim(x_min - 0.5, x_max + 0.5)
                    ax.set_ylim(y_min - 0.5, y_max + 0.5)
                    ax.grid(False)
                    ax.set_title(fr"${t}$-plane")
                    ax.set_xlabel(r"${x}$")
                    ax.set_ylabel(r"${y}$")

                axs[0,3].set_xlabel(r"${x}$")
                axs[0,3].set_ylabel(r"${y}$")

                fig.suptitle(rf"{run_number}   $\beta = {beta}$")

                #print(np.mean(abs(config[:,2])))

                n = int(np.sqrt(len(config[:, 0])))
                assert n**2 == len(config[:, 0]), "The length of config[:,0] must be a perfect square."

                # Create q-values directly over the [-pi, pi) range
                qx = 2*np.pi/n*np.arange(n)
                qy = 2*np.pi/n*np.arange(n)

                points = np.array([[x, y, 0] for y in qy for x in qx])

                b1 = 2*np.pi*np.array([1, 0, 0])
                b2 = 2*np.pi*np.array([0, 1, 0])
                N = len(points[:,0])
                m = 4
                qs = np.tile(points, [m,1]) + np.concatenate((np.zeros(np.shape(points)), np.tile(-b1, [N,1]), np.tile((-b1-b2), [N,1]), np.tile(-b2, [N,1])))

                QX = qs[:,0]
                QY = qs[:,1]
                #QX, QY = np.meshgrid(qx, qy)

                Sx = scipy.fft.fft2(config[:, 0].reshape(n, n))
                Sy = scipy.fft.fft2(config[:, 1].reshape(n, n))
                Sz = scipy.fft.fft2(config[:, 2].reshape(n, n))

                SQ = np.real(np.conj(Sx) * Sx + np.conj(Sy) * Sy + np.conj(Sz) * Sz).reshape(-1)

                SQold = SQ
                for i in range(m-1):
                    SQ = np.concatenate((SQ, SQold))

                plt.gca().set_aspect('equal')
                triang = tri.Triangulation(QX, QY)
                sqplot = axs[1,0].tricontourf(triang, SQ, 100, alpha = 1)
                divider = make_axes_locatable(axs[1,0])
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(sqplot, cax=cax)
                axs[1,0].plot([-np.pi, -np.pi, np.pi, np.pi, -np.pi], [-np.pi, np.pi, np.pi, -np.pi, -np.pi], "k")
                lims = [-0.9, 0.9]
                axs[1,0].set_xlim(lims)
                axs[1,0].set_ylim(lims)
                axs[1,0].plot(Qs1[:,0], Qs1[:,1], "--", color='k', linewidth=0.5, alpha=0.5)
                axs[1,0].set_xlabel(r"$q_x$")
                axs[1,0].set_ylabel(r"$q_y$")

                confQ, confn, confplanarity = BasisFunctions.ComputeQ_n_planarity(config, nbrs)

                confQtheta = np.mod(confQ[:,1],np.pi)
                confQr = confQ[:,0]

                print(beta)
                print(np.min(confQtheta),np.max(confQtheta))

                colorsconfQ = colors(np.cos(2*confQtheta), np.sin(2*confQtheta), np.zeros(len(confQ[:,0])))
                colorsconfn = colors(confn[:,0],confn[:,1],confn[:,2])
                colorsconfQr = (confQr)

                axs[1,1].imshow(colorsconfQ , aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))

                axs[1,2].imshow(colorsconfn , aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))

                Qrplot = axs[1,3].imshow(np.flip(colorsconfQr.reshape((len(y_unique), len(x_unique))), axis=0), cmap='hot', aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))
                #Qrplot = axs[1,3].imshow(np.flip(confn[:,2].reshape((len(y_unique), len(x_unique))), axis=0), cmap='hot', aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))
                divider = make_axes_locatable(axs[1,3])
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(Qrplot, cax=cax)

                planarityplot = axs[0,3].imshow(np.flip(confplanarity.reshape((len(y_unique), len(x_unique))), axis=0), cmap='hot', aspect='equal', extent=(x_min - 0.5, x_max + 0.5, y_min - 0.5, y_max + 0.5))
                divider = make_axes_locatable(axs[0,3])
                cax = divider.append_axes("right", size="5%", pad=0.05)
                fig.colorbar(planarityplot, cax=cax)

                for ax, t in zip(axs[1,:3], ["S(\vec{q})", "\hat{q}", "\hat{n}"]):
                    ax.grid(False)

            
                axs[0,3].set_title(r"$\vec{S}_{\vec{r}} \cdot (\vec{S}_{\vec{r}+\hat{x}} \times \vec{S}_{\vec{r}+\hat{y}})$")
                
                axs[1,0].set_title(r"$S(\vec{q})$")
                axs[1,1].set_title(r"$\hat{q}$")
                axs[1,2].set_title(r"$\hat{n}$")
                axs[1,3].set_title(r"$|\vec{q}|$")

                for ax in axs[1,1:4]:
                    ax.set_xlabel(r"${x}$")
                    ax.set_ylabel(r"${y}$")

                    ax.set_xlim(x_min - 0.5, x_max + 0.5)
                    ax.set_ylim(y_min - 0.5, y_max + 0.5)

                minz = np.sort(confn[:,2])
                print(minz[:100])

                """fig2, ax2 = plt.subplots(1, 2, figsize=(10, 6))

                ax2.plot()


                fig2.suptitle(rf"{run_number}   $\beta = {beta}$")"""
                

    plt.show()