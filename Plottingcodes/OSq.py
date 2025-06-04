import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

import BasisFunctions

import plotstyle

try:
    a = eval(sys.argv[1])
except:
    a = eval(input('Multiplication factor: '))

print(a)

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

startb = 1
stopb  = 10
bins = range(startb,stopb+1)

ibdata = get_data("indxbeta.txt", ["runs", "betas"])

runs = ibdata["runs"]
betas = ibdata["betas"]

Nb = 0
Nq = 64
#Nq = 128
"""for line in infile:
    words = line.split()
    elif len(words) == 15 and Nb == 1:
        Nq += 1

infile.close()"""

print(Nq)

Sqs = np.zeros((Nq, len(betas)))

for bin in bins:
    ofs = ['qcorr_' + str(run) + '_' + str(bin) + '.dat' for run in runs]
    for (file,ib) in zip(ofs, range(len(betas))):
        qs = []
        Sq = []
        infile = open(file, 'r')
        for line in infile:
            words = [float(f) for f in line.split()]
            qs.append(np.array(words[:3]))
            Sq.append(words[-1])

        Sqs[:, ib] += Sq

Sqs = np.array(Sqs)/len(bins)
qs = np.array(qs)

folder0 = os.path.expanduser('~/Documents/MCHeisenberg/Data/')
T0files = ['Pyrochlore32/NoHole_MINE/Run007']
labels = ['cyc=1,bins=1e7,dE=1e-10']

T0files = []
labels0 = []

"""
if(lattice == "PYROCHLORE"):
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
"""

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

qinds = [0,1,2,3,4,10,11,12,13,20,21,34,42]
qinds = range(Nq)

fig, ax = plt.subplots()

print(1./np.array(betas)[1:3])
#for i in range(int(Nq/2+1)):
for i in qinds: #range(len(qs)):
    p = np.polyfit(1./np.array(betas)[1:-1], a*Sqs[i, :][1:-1], 2)
    Temps = np.linspace(0, 1, 100)
    color=next(ax._get_lines.prop_cycler)['color']
    #color = 'b'
    plt.semilogx(np.array(betas), a*Sqs[i, :], 'o-', markersize = 3, c = color)#, label = labels[i])
    #plt.plot(Temps, p[0]*Temps**2 + p[1]*Temps + p[2], '-', markersize = 3, c = color, alpha = 0.5)#, label = labels[i])

    #plt.plot(1./np.array(betas), Sqs[i, :], 'o-', markersize = 3, c = color)#, label = labels[i])
    for j in range(len(T0files)):
        plt.plot([0, 1], [SqT0[j, i],SqT0[j, i]], '--', c = color, label=labels[j])
    plt.xlabel(r"$\beta$ ", fontsize=14)
    plt.ylabel(r"$S^{zz}(Q)$", fontsize=14)
    #plt.xlim([0,1000])
    #plt.ylim(ylims)
    #plt.legend()

BasisFunctions.plotter(qs, Sqs[:,-1],"hhl","Pyrochlore")
BasisFunctions.plotter(qs, Sqs[:,-1],"hk0","Pyrochlore")

for j in range(len(T0files)):
    BasisFunctions.plotter(qs, SqT0[j, :],"hhl","Pyrochlore")
    BasisFunctions.plotter(qs, SqT0[j, :],"hk0","Pyrochlore")

plt.show()


sum = np.zeros(len(betas))
for i in range(len(betas)):
    sum[i] = np.sum(Sqs[:,i])

print(np.array(betas))
