import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os
import sys

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

Lx = int(paramfile.readline().split()[-1])
Ly = int(paramfile.readline().split()[-1])
Lz = int(paramfile.readline().split()[-1])

a = " "
while a != 'nbins':
    words = paramfile.readline().split()
    try:
        a = words[0]
        b = words[-1]
    except:
        continue

nbins = int(b)
nbins = 100

print(a, nbins)

print(Lx, Ly, Lz)

if lattice == "PYRO16":
    Nq = 64*Lx*Ly*Lz

qs = []
Sqmean = []
infile = open("Run" + run_number + '/expSq.txt', 'r')
infile.readline()
for line in infile:
    words = [float(f) for f in line.split()]
    qs.append(np.array(words[:3]))
    Sqmean.append(words[3])
    print(words[3])

print(len(Sqmean))

qs = np.array(qs)
Sqmean = np.array(Sqmean)

infile = open("Run" + run_number + "/expSqlog.txt", 'r')

SqT0 = np.zeros((nbins,Nq))

for i in range(nbins):
    for j in range(Nq):
        infile.readline()
        infile.readline()
        SqT0[i,j] = float(infile.readline().split()[3])

Sqmin = np.min(SqT0)
Sqmax = np.max(SqT0)

plt.figure()
plt.hist(SqT0[0:100, 0])
plt.xlim([Sqmin,Sqmax])
plt.show()

"""
for j in range(8):
    fig, ax = plt.subplots(Lx*Ly*Lz,8)
    #for i in range(int(Nq/2+1)):
    for i in range(8*Lx*Ly*Lz):
        ax[i//8,i%8].hist(SqT0[:,j*8*Lx*Ly*Lz + i])
        plt.xlabel(r"$S^{zz}(Q)$", fontsize=14)
        plt.ylabel(r"Count", fontsize=14)
        plt.xlim([Sqmin,Sqmax])
        #plt.ylim(ylims)
        #plt.legend()
"""

plt.show()
