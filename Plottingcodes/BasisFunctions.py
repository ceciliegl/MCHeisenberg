import numpy as np
import scipy
import scipy.interpolate
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt

import sys

import matplotlib.tri as tri
import cmath

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

mpl.rc('xtick', labelsize=14)
mpl.rc('ytick', labelsize=14)
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

mpl.rcParams.update({'font.size': 14})

def SzqSzq(q0, ZZ, r, totsites, p):
    s = 0
    for i in range(totsites):
        for j in range(totsites):
            a = ZZ[i][j]*np.exp(-1j*np.dot(q0, r[i]-r[j])) #Should be possible to vectorize this?
            s += a

    return s/(totsites)



def plotter(qq, MM, plottype, latticetype):
    if latticetype == "KAGOME":
        plottype = "hk0"

    if plottype == "hhl":
        plotqs = []
        plotMags = []

        for k in range(len(qq)):
            if abs(qq[k][0]-qq[k][1]) < 1e-6:
                plotqs.append(qq[k])
                plotMags.append(MM[k])

    elif plottype == "hk0":
        plotqs = []
        plotMags = []

        for k in range(len(qq)):
            if abs(qq[k][2]) < 1e-6:
                plotqs.append(qq[k])
                plotMags.append(MM[k])

    elif plottype == "hkk":
        plotqs = []
        plotMags = []

        for k in range(len(qq)):
            if abs(qq[k][1]-qq[k][2]) < 1e-6:
                plotqs.append(qq[k])
                plotMags.append(MM[k])

    elif plottype == "0kl":
        plotqs = []
        plotMags = []

        for k in range(len(qq)):
            if abs(qq[k][0]) < 1e-6:
                plotqs.append(qq[k])
                plotMags.append(MM[k])

    plotqs = np.array(plotqs)

    Z = np.array(plotMags)
    if latticetype == "KAGOME":
        b1 = 2*2*np.pi*np.array([1,-1/np.sqrt(3),0])
        b2 = 2*2*np.pi*np.array([0,2./np.sqrt(3),0])
        N = len(plotqs[:,0])
        n = 14
        qs = np.tile(plotqs, [n,1]) + np.concatenate((np.zeros(np.shape(plotqs)), np.tile(-b1, [N,1]), np.tile((-b1-b2), [N,1]), np.tile(-b2, [N,1]), np.tile(b1, [N,1]), np.tile((b1+b2), [N,1]), np.tile(b2, [N,1]), np.tile(-2*b1, [N,1]), np.tile(-2*b2, [N,1]), np.tile(-2*b1+b2, [N,1]), np.tile(-2*b2+b1, [N,1]), np.tile(-2*b1-b2, [N,1]), np.tile(-2*b2-b1, [N,1]), np.tile(-2*(b1+b2), [N,1])))
    elif latticetype == "PYRO16":
        if plottype == "hhl":
            b1 = 8*np.pi*np.array([1, 1, 0])
            b2 = 8*np.pi*np.array([0, 0, 1])
        if plottype == "hk0":
            b1 = 8*np.pi*np.array([1, 0, 0])
            b2 = 8*np.pi*np.array([0, 1, 0])
        if plottype == "hkk":
            b1 = 8*np.pi*np.array([1, 0, 0])
            b2 = 8*np.pi*np.array([0, 1, 1])
        if plottype == "0kl":
            b1 = 8*np.pi*np.array([0, 1, 0])
            b2 = 8*np.pi*np.array([0, 0, 1])
        N = len(plotqs[:,0])
        n = 9
        qs = np.tile(plotqs, [n,1]) + np.concatenate((np.zeros(np.shape(plotqs)), np.tile(-b1, [N,1]), np.tile((-b1-b2), [N,1]), np.tile(-b2, [N,1]), np.tile(b1, [N,1]), np.tile(b2, [N,1]), np.tile((b1+b2), [N,1]), np.tile((-b1+b2), [N,1]), np.tile((b1-b2), [N,1])))
    else:
        n = 1
        qs = plotqs

    Zold = Z
    for i in range(n-1):
        Z = np.concatenate((Z, Zold))


    if plottype == "hhl":
        X = qs[:,0]
        Y = qs[:,2]
    elif plottype == "hk0":
        X = qs[:,0]
        Y = qs[:,1]
    elif plottype == "hkk":
        X = qs[:,0]
        Y = qs[:,2]
    elif plottype == "0kl":
        X = qs[:,1]
        Y = qs[:,2]

    fig = plt.figure()
    plt.gca().set_aspect('equal')

    """a = int(np.sqrt(len(X)))

    xx = np.zeros((a, a))
    yy = np.zeros((a, a))
    zz = np.zeros((a, a))

    for i in range(len(X)):
        print(i//a, i%a)
        xx[i//a,i%a] = X[i]
        yy[i//a,i%a] = Y[i]
        zz[i//a,i%a] = Z[i]

    #print(xx)

    plt.scatter(X, Y, s=200, c=Z)"""

    if latticetype == 'KAGOME':
        triang = tri.Triangulation(X, Y)
        plt.tricontourf(triang, Z, 100, alpha = 1)
    else:
        xx,yy = np.meshgrid(np.unique(X),np.unique(Y))
        rbf = scipy.interpolate.Rbf(X, Y, Z, function='linear')
        zz = rbf(xx, yy)

        plt.contourf(xx, yy, zz, 100, alpha = 1)

    if latticetype == 'KAGOME':
        Kx = np.array([4*np.pi/3, 2*np.pi/3, -2*np.pi/3, -4*np.pi/3, -2*np.pi/3, 2*np.pi/3, 4*np.pi/3])
        Ky = np.array([0, -2*np.pi/np.sqrt(3), -2*np.pi/np.sqrt(3), 0, 2*np.pi/np.sqrt(3), 2*np.pi/np.sqrt(3), 0])
        plt.plot(Kx, Ky, color = "k", linewidth = 0.5)
        plt.plot(Kx+b1[0]/2*np.ones(7), Ky+b1[1]/2*np.ones(7), color = "k", linewidth = 0.5)
        plt.plot(Kx+(b1+b2)[0]/2*np.ones(7), Ky+(b1+b2)[1]/2*np.ones(7), color = "k", linewidth = 0.5)
        plt.plot(Kx+b2[0]/2*np.ones(7), Ky+b2[1]/2*np.ones(7), color = "k", linewidth = 0.5)
        plt.plot(Kx+(-b1)[0]/2*np.ones(7), Ky+(-b1)[1]/2*np.ones(7), color = "k", linewidth = 0.5)
        plt.plot(Kx+(-b1-b2)[0]/2*np.ones(7), Ky+(-b1-b2)[1]/2*np.ones(7), color = "k", linewidth = 0.5)
        plt.plot(Kx+(-b2)[0]/2*np.ones(7), Ky+(-b2)[1]/2*np.ones(7), color = "k", linewidth = 0.5)

        c = 4*np.pi
        plt.xlim([-c,c])
        plt.ylim([-c,c])
        plt.xticks([-4*np.pi, -3*np.pi, -2*np.pi, -np.pi, 0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi], [r"$-4\pi$", r"$-3\pi$", r"$-2\pi$", r"$-\pi$", r"$0$", r"$\pi$", r"$2\pi$", r"$3\pi$", r"$4\pi$"])
        plt.yticks([-4*np.pi, -3*np.pi, -2*np.pi, -np.pi, 0, np.pi, 2*np.pi, 3*np.pi, 4*np.pi], [r"$-4\pi$", r"$-3\pi$", r"$-2\pi$", r"$-\pi$", r"$0$", r"$\pi$", r"$2\pi$", r"$3\pi$", r"$4\pi$"])
    elif latticetype == 'PYRO16':
        c = 8*np.pi
        plt.xlim([-c,c])
        plt.ylim([-c,c])
        plt.xticks([-8*np.pi, -4*np.pi, 0, 4*np.pi, 8*np.pi], [r"$-8\pi$", r"$-4\pi$", r"$0$", r"$4\pi$", r"$8\pi$"])
        plt.yticks([-8*np.pi, -4*np.pi, 0, 4*np.pi, 8*np.pi], [r"$-8\pi$", r"$-4\pi$", r"$0$", r"$4\pi$", r"$8\pi$"])

    if plottype == "hhl":
        plt.title(r'$(hhl)$')
        plt.xlabel(r'$h$')
        plt.ylabel(r'$l$')

    elif plottype == "hk0":
        plt.title(r'$(hk0)$')
        plt.xlabel(r'$h$')
        plt.ylabel(r'$k$')

    elif plottype == "hkk":
        plt.title(r'$(hkk)$')
        plt.xlabel(r'$h$')
        plt.ylabel(r'$k$')

    elif plottype == "0kl":
        plt.title(r'$(0kl)$')
        plt.xlabel(r'$k$')
        plt.ylabel(r'$l$')

    plt.colorbar()

def GenLattice(latticetype, Lx, Ly, Lz):
    if latticetype == 'CHAIN':
        a1 = np.array([1,0,0])
        a2 = np.array([0,1,0])
        a3 = np.array([0,0,1])

        nsubl = 1
        alphas = np.array([[0,0,0]])

        dim = 1
    elif latticetype == 'TRIANGULAR':
        a1 = np.array([0,1,0])
        a2 = np.array([np.sqrt(3)/2,0.5,0])
        a3 = np.array([0,0,1])

        nsubl = 1
        alphas = np.array([[0,0,0]])

        dim = 2
    elif latticetype == 'KAGOME':
        a1 = np.array([0,1,0])
        a2 = np.array([np.sqrt(3)/2,0.5,0])
        a3 = np.array([0,0,1])

        nsubl = 3
        alphas = np.array([[0,0,0],0.5*a1,0.5*a2])

        dim = 2
    elif latticetype == 'PYROCHLORE':
        a1 = np.array([0,1/2.,1/2.])
        a2 = np.array([1/2.,0,1/2.])
        a3 = np.array([1/2.,1/2.,0])

        nsubl = 4
        alphas = np.array([[0,0,0],0.5*a1,0.5*a2,0.5*a3])

        dim = 3
    elif latticetype == 'PYRO16':
        a1 = np.array([1,0,0])
        a2 = np.array([0,1,0])
        a3 = np.array([0,0,1])
        nsubl = 16
        dim = 3

        c1 = (a2+a3)/2; c2 = (a3+a1)/2; c3 = (a1+a2)/2;

        alphas = np.array([[0,0,0],c1/2,c2/2,c3/2, c1, c1+c1/2, c1+c2/2, c1+c3/2, c2, c2+c1/2, c2+c2/2, c2+c3/2, c3, c3+c1/2, c3+c2/2, c3+c3/2])

    else:
        print('WARNING: a vectors not defined for latticetype ', latticetype)

    Ltot = Lx*Ly*Lz
    totsites = Ltot*nsubl

    b1 = 2*np.pi*np.cross(a2,a3)/np.dot(a1, np.cross(a2,a3));
    b2 = 2*np.pi*np.cross(a3,a1)/np.dot(a2, np.cross(a3,a1));
    b3 = 2*np.pi*np.cross(a1,a2)/np.dot(a3, np.cross(a1,a2));

    r = np.zeros((totsites,3))
    j = 0
    for iz in range(Lz):
        for iy in range(Ly):
            for ix in range(Lx):
                for isubl in range(nsubl):
                    r[j] = ix*a1 + iy*a2 + iz*a3 + alphas[isubl]
                    j += 1

    if latticetype == 'KAGOME':
        nZones = 4
        zones=[np.array([0,0,0]),b1, b2, b1+b2]
    elif latticetype == 'PYROCHLORE':
        nZones = 8
        zones=[np.array([0,0,0]),b1, b2, b3, b1+b2, b2+b3, b3+b1, b1+b2+b3]
    elif latticetype == 'PYRO16':
        nZones = 4*4*4
        zones=[]
        for i1 in range(4):
            for i2 in range(4):
                for i3 in range(4):
                    zones.append(i1*b1+i2*b2+i3*b3)
    else:
        zones = [np.array([0,0,0])]
        nZones = 1

    q = np.zeros((Ltot*nZones,3))
    j = 0

    #print(zones)

    for zone in zones:
        for iz in range(Lz):
            for iy in range(Ly):
                for ix in range(Lx):
                    q[j] = ix/float(Lx)*b1 + iy/float(Ly)*b2 + iz/float(Lz)*b3 + zone
                    j += 1

    return r, q, totsites

def ComputeQ_n_planarity(config, nbrs):
    ansQ = np.zeros((len(config),2)) #We save this in r, phi format.
    ansn = np.zeros((len(config), 3)) #We save this as nx, ny, nz.

    planarity = np.zeros(len(config))

    n0 = np.array([0,0,1])

    for i in range(len(config)):
        spin0 = config[i]
        spin1 = config[nbrs[i][0]]
        spin2 = config[nbrs[i][1]]

        cosQa1 = np.dot(spin0,spin1)
        cosQa2 = np.dot(spin0,spin2)

        cross01 = np.cross(spin0,spin1)
        cross02 = np.cross(spin0,spin2)


        n = cross01/np.linalg.norm(cross01)
        m = cross02/np.linalg.norm(cross02) #Is this very different from n?
        sinQa1 = np.linalg.norm(cross01)
        sinQa2 = np.linalg.norm(cross02)

        if np.dot(n, m) < 0:
            m = -m
            sinQa2 = -sinQa2

        """
        if n[2] < 0:
            n = -n
            sinQa1 = -sinQa1

        if m[2] < 0:
            m = -m
            sinQa2 = -sinQa2
        """
        
        planarity[i] = np.abs(np.dot(spin0,np.cross(spin1,spin2)))


        Qa1 = np.arctan2(sinQa1, cosQa1)
        Qa2 = np.arctan2(sinQa2, cosQa2)

        Q = np.array([Qa1, Qa2]) #ONLY FOR SQUARE LATTICE
        ansQ[i] = [np.linalg.norm(Q), np.arctan2(Q[1],Q[0])]

        #print(np.arctan2(Q[1],Q[0]))

        ansn[i] = n

        n0 = n

    return ansQ, ansn, planarity

def coord_to_site(coord, Nx, Ny, Nz=1):
  #Store sites as (Nx*(Nz*(nz)+ny)+nx)
  nx = (coord[0]+Nx)%Nx
  ny = (coord[1]+Ny)%Ny
  nz = (coord[2]+Nz)%Nz
  return int(Nx*(Ny*(nz)+ny)+nx)