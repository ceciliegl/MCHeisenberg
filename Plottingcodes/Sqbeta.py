import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

import BasisFunctions

import plotstyle

try:
    run_number = sys.argv[1]
except:
    run_number = input('Run number: ')

import matplotlib
matplotlib.rc('xtick', labelsize=14)
matplotlib.rc('ytick', labelsize=14)
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def get_data(filename, variables):
    df = pd.read_csv(filename,\
                    delim_whitespace=True, \
                    engine='python', \
                    names=variables)
    return df
    #using pandas to read the data files

if float(run_number) < 10:
    run_number = "00" + run_number
elif float(run_number) < 100:
    run_number = "0" + run_number
elif float(run_number) < 1000:
    run_number = run_number
else:
    print("Run number too big")
    exit(1)

paramfile = open("Run" + run_number + "/parameters.txt", 'r')

lattice = paramfile.readline().split()[-1]

print(lattice)

Nx = int(paramfile.readline().split()[-1])
Ny = int(paramfile.readline().split()[-1])
Nz = int(paramfile.readline().split()[-1])


infile = open("Run" + run_number + "/expSq.txt", 'r')

Nb = 0
Nq = 0
for line in infile:
    words = line.split()
    if len(words) == 1:
        Nb += 1
    elif len(words) == 15 and Nb == 1:
        Nq += 1

infile.close()

print(Nq)

Szzq = np.zeros((Nq, Nb))
dSzzq = np.zeros((Nq, Nb))
betas = np.zeros(Nb)
qs = np.zeros((Nq, 3))

infile = open("Run" + run_number + "/expSq.txt", 'r')

Nb = 0
Nq = 0
for line in infile:
    words = line.split()
    if len(words) == 1:
        Nb += 1
        Nq = 0
        betas[Nb-1] = float(words[0])
    elif len(words) == 15:
        Nq += 1
        Szzq[Nq-1, Nb-1] = float(words[-4])
        dSzzq[Nq-1, Nb-1] = float(words[-3])
        qs[Nq-1, 0] = float(words[0])
        qs[Nq-1, 1] = float(words[1])
        qs[Nq-1, 2] = float(words[2])

print(betas)
print(qs)
print(Nq)

if(lattice == "KAGOME"):
    qinds = range(Nq)
    labels = ["Q = (%.2f, %.2f, %.2f)" % (qs[i, 0], qs[i, 1], qs[i, 2]) for i in qinds]
    ylims = [0,1]
if(lattice == "PYROCHLORE"):
    qinds = range(Nq)
    labels = ["Q = (%.2f, %.2f, %.2f)" % (qs[i, 0], qs[i, 1], qs[i, 2]) for i in qinds]
    ylims = [0,1]
if(lattice == "PYRO16"):
    qinds = range(Nq)
    labels = ["Q = (%.2f, %.2f, %.2f)" % (qs[i, 0], qs[i, 1], qs[i, 2]) for i in qinds]
    #qinds = [0,1,4,6,16,19,28]
    #labels = [r"$\Gamma$",r"$\frac{1}{2}b_i$",r"$b_i$",r"$\frac{1}{2}b_i+b_j$",r"$b_i + b_j$",r"$\frac{1}{2}b_i + b_j + b_k$",r"$b_1 + b_2 + b_3$"]
    ylims = [0,1]
elif(lattice == "TRIANGULAR"):
    qinds = range(Nq)
    labels = ["Q = (%.2f, %.2f, %.2f)" % (qs[i, 0], qs[i, 1], qs[i, 2]) for i in qinds]
    ylims = [0,3]
    #myqs = [0, 1, 5]
    #labels = [r"$\Gamma$", r"$q$", r"K"]
else:
    qinds = range(int(Nq/2+1))
    labels = ["Q = (%.2f, %.2f, %.2f)" % (qs[i, 0], qs[i, 1], qs[i, 2]) for i in qinds]


folder0 = os.path.expanduser('~/Documents/MCHeisenberg/Data/')

if(lattice == "KAGOME"):
    if Nx == 2 and Ny ==2 and Nz == 1:
        T0files = ['Kagome12/NoHole_MINE/Run006']
        labels0 = ['cyc=1,bins=1e6,dE=1e-10']
    elif Nx == 3 and Ny ==3 and Nz == 1:
        T0files = ['Kagome27/NoHole_MINE/Run006']
        labels0 = ['cyc=1,bins=1e6,dE=1e-10']
    else:
        T0files = []
        labels0 = []
elif(lattice == "PYROCHLORE"):
    if Nx == 2 and Ny ==2 and Nz == 2:
        T0files = ['Pyrochlore32/NoHole_MINE/Run007']
        labels0 = ['cyc=1,bins=1e7,dE=1e-10']
    else:
        T0files = []
        labels0 = []
elif(lattice == "PYRO16"):
    if Nx == 1 and Ny == 1 and Nz == 1:
        T0files = ['PYRO16_MINE_Varynbins/NoHole/Run005']
        labels0 = ['cyc=1,bins=1e5,dE=1e-10']
    elif Nx == 2 and Ny ==2 and Nz == 1:
        T0files = ['PYRO64/NoHole_MINE_Varynbins/Run006']
        labels0 = ['cyc=1,bins=1e6,dE=1e-10']
    elif Nx == 2 and Ny ==2 and Nz == 2:
        T0files = ['PYRO128/NoHole_MINE_Varynbins/Run006']
        labels0 = ['cyc=1,bins=1e6,dE=1e-10']
    elif Nx == 8 and Ny ==1 and Nz == 1:
        T0files = ['PYRO128/NoHoleASym_MINE_Varynbins/Run005']
        labels0 = ['cyc=1,bins=1e5,dE=1e-10']
    else:
        T0files = []
        labels0 = []
else:
    T0files = []
    labels0 = []

SqT0 = np.zeros((len(T0files),Nq))
for i, T0file in enumerate(T0files):
    qs = []
    Sq = []
    infile = open(folder0 + T0file + '/expSq.txt', 'r')
    infile.readline()
    for line in infile:
        words = [float(f) for f in line.split()]
        qs.append(np.array(words[:3]))
        Sq.append(words[3])
        print(words[3])

    print(len(Sq))

    SqT0[i,:] += np.array(Sq)


"""olavruns = [0,7,10,15,17,22,25,30,32,37,40,45]
startb = 1
stopb  = 10
bins = range(startb,stopb+1)
olavbetas = [1e-1,1,2,1e1,2e1,1e2,2e2,1e3,2e3,1e4,2e4,1e5]

olavSqs = np.zeros((Nq, len(olavbetas)))

for bin in bins:
    ofs = ['qcorr_' + str(olavrun) + '_' + str(bin) + '.dat' for olavrun in olavruns]
    for (file,ib) in zip(ofs, range(len(olavbetas))):
        olavqs = []
        Sq = []
        infile = open(file, 'r')
        for line in infile:
            words = [float(f) for f in line.split()]
            olavqs.append(np.array(words[:3]))
            Sq.append(words[-1])

        olavSqs[:, ib] += Sq

olavSqs = np.array(olavSqs)/len(bins)
olavqs = np.array(olavqs)"""

"""T0files = ['qcorr_T0_1.dat', 'qcorr_T0_2.dat']

SqT0 = np.zeros(Nq)
for T0file in T0files:
    olavqs = []
    Sq = []
    infile = open(T0file, 'r')
    for line in infile:
        words = [float(f) for f in line.split()]
        olavqs.append(np.array(words[:3]))
        Sq.append(words[-1])
        print(words[-1])

    print(len(Sq))

    SqT0 += np.array(Sq)

SqT0 /= len(T0files)"""

fig, ax = plt.subplots()
#for i in range(int(Nq/2+1)):
"""
for i in range(len(olavqs)):
    color=next(ax._get_lines.prop_cycler)['color']
    color = 'b'
    plt.semilogx(np.array(olavbetas), olavSqs[i, :], 'o-', markersize = 3, c = color)#, label = labels[i])
    #plt.plot(1./np.array(olavbetas), olavSqs[i, :], 'o-', markersize = 3, c = color)#, label = labels[i])
    #plt.plot([0], [SqT0[i]], 'bo')
    plt.xlabel(r"$\beta$ ", fontsize=14)
    plt.ylabel(r"$S^{zz}(Q)$", fontsize=14)
    #plt.xlim([0,1000])
    #plt.ylim(ylims)
    #plt.legend()
"""

#fig, ax = plt.subplots()
#for i in range(int(Nq/2+1)):

for i in range(len(qinds)):
    color=next(ax._get_lines.prop_cycler)['color']
    #color = 'magenta'
    plt.semilogx(betas, Szzq[qinds[i],:], 'o-', markersize = 3, c = color, label = labels[i])
    plt.errorbar(betas, Szzq[qinds[i],:], dSzzq[qinds[i],:], capsize=3, elinewidth=1.5, markeredgewidth=1.5, c = color)
    #plt.plot(1./betas, Szzq[qinds[i],:], 'o-', markersize = 3, c = color, label = labels[i])
    plt.xlabel(r"$T$ ", fontsize=14)
    plt.ylabel(r"$S^{zz}(Q)$", fontsize=14)
    for j in range(len(T0files)):
        plt.plot([0, 100000], [SqT0[j, i],SqT0[j, i]], '--', c = color, label=labels0[j])
    #plt.xlim([0,1000])
    #plt.ylim(ylims)
    #plt.legend()

BasisFunctions.plotter(qs, Szzq[:,-1],"hhl",lattice)
BasisFunctions.plotter(qs, Szzq[:,-1],"hk0",lattice)

for j in range(len(T0files)):
    BasisFunctions.plotter(qs, SqT0[j, :],"hhl",lattice)
    BasisFunctions.plotter(qs, SqT0[j, :],"hk0",lattice)

plt.show()
