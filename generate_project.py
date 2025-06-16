import os
import numpy as np
import random

mainproject = "SSL_Square_RING200x200"  #Set to zero if only one project.
project = "SI0_minE"
description = "Running SSL square lattice."
jobname = "myjob"
time = "5:00:00"
runmin = 10
runmax = 10
runsame = 0
nruns = (runmax-runmin) + 1
NICE = 11

#LIBS#
#BOOST = 0       #Higher precision in Eigen-calculations. Time-consuming. Not implemented for now.

#LATTICE#
lattice   = "SQUARE"
NX        = 200*np.ones(nruns, int) #np.array([2], int) #2*np.ones(nruns, int)
nruns     = len(NX)
runmax    = runmin + (nruns-1)
NY        = 200*np.ones(nruns, int)
if lattice == "FCC" or lattice == "DIAMOND" or lattice == "TETRAGONALDIAMOND" or lattice == "PYROCHLORE": #If lattice is three dimensional...
    NZ    = 1*np.ones(nruns, int)
else:
    NZ    = np.ones(nruns, int)

Nh = 0*np.ones(nruns, int)

holepos = runmin

OBC = 0

#EXCHANGE#
J1     = -1*np.ones(nruns)
J2     = 0.28*np.ones(nruns)
J3     = J2/2. #0.25*np.ones(nruns)#np.array([-1.,-1./2,-1./3,-1./5,-1./8,-1./10,-1./15,-1./20,-1./30,0,1./30,1./20,1./15,1./10,1./8,1./5,1./3,1./2,1.])
J3a    = 0*np.ones(nruns)
J3b    = 0*np.ones(nruns)
DXY    = 0*np.ones(nruns) #Single ion anisotropy, favouring XY-plane.

ncyc        = int(1e1)
nequ        = int(1e3)
nbins       = int(1)
nconew      = int(100)
print       = int(ncyc)    #Write out data every print cycle.
configprint = int(1)       #Write out spin configuration at end of each beta.


ncyc        = int(ncyc)
nequ        = int(nequ)
nbins       = int(nbins)
nconew      = int(nconew)
print       = int(print)
configprint = int(configprint)

RANSEED = np.random.randint(1e8,1e10,nruns)

EQUENERGY = 123456789 #If minimize energy, measuring will start when this energy is reached. Set to 123456789 if you just want to print out the spin configuration after running the program. (123456789 NOT IMPLEMENTED)

RESETOLDFILES = 1



if mainproject:
    os.system("mkdir " + "Data/" + mainproject)
    totalproject = mainproject + "/" + project
else:
    totalproject = project

os.system("mkdir " + "Data/" + totalproject)

descfile = open("Data/" + totalproject + "/description.txt", 'w')
descfile.write(description)

runsub = open("Data/" + totalproject + "/run.sub", "w")

runsub.write("#!/bin/sh\n")
runsub.write("#SBATCH --account=nn4563k\n")
runsub.write("#SBATCH --time=" + time + "\n")
runsub.write("#SBATCH --partition=bigmem\n")
runsub.write("#SBATCH --mem-per-cpu=16999M\n")
runsub.write("#SBATCH --ntasks=1\n")
runsub.write("#SBATCH --signal=B:USR1@60\n")

if runsame:
    runsub.write("#SBATCH --array={}-{}%{}\n".format(runmin, runmax, runsame))
else:
    runsub.write("#SBATCH --array={}-{}\n".format(runmin, runmax))
runsub.write("#SBATCH --job-name="+jobname+"\n")
runsub.write("#SBATCH --error=stderr.dat\n")
runsub.write("#SBATCH --output=stdout.dat\n")
runsub.write("#set -o errexit\n")
runsub.write("\n")
runsub.write("\n")
runsub.write('#cleanup "rsync -av $SCRATCH/ $SUBMITDIR/ --exclude=stdout.dat --exclude=stderr.dat"\n')
runsub.write("#cd $SUBMITDIR\n")
runsub.write("#rsync -av $SUBMITDIR/ $SCRATCH/ --exclude=rundir\n")
runsub.write("#cd $SCRATCH\n")
runsub.write("echo Running program.....\n")
runsub.write("$HOME/Documents/MCHeisenberg/Code/program " + totalproject + " $SLURM_ARRAY_TASK_ID\n")





for run in range(runmin, nruns + runmin):
    #delta = np.ones(num)*(np.logspace(dmin, dmax, num)[run])
    if run < 10:
        os.system("mkdir " + "Data/" + totalproject + "/Run00" + str(run))
        outfile = open("Data/" + totalproject + "/Run00" + str(run) + "/parameters.txt", 'w')
    elif run < 100:
        os.system("mkdir " + "Data/" + totalproject + "/Run0" + str(run))
        outfile = open("Data/" + totalproject + "/Run0" + str(run) + "/parameters.txt", 'w')
    elif run < 1000:
        os.system("mkdir " + "Data/" + totalproject + "/Run" + str(run))
        outfile = open("Data/" + totalproject + "/Run" + str(run) + "/parameters.txt", 'w')
    else:
        print("generate project: Run number is bigger than 999")
        exit(1)

    outfile.write("Lattice = ")
    outfile.write(lattice)
    outfile.write("\n")
    outfile.write("NX = ")
    outfile.write(str(NX[run-runmin]))
    outfile.write("\n")
    outfile.write("NY = ")
    outfile.write(str(NY[run-runmin]))
    outfile.write("\n")
    outfile.write("NZ = ")
    outfile.write(str(NZ[run-runmin]))
    outfile.write("\n")
    outfile.write("Nh = ")
    outfile.write(str(Nh[run-runmin]))
    outfile.write("\n")
    outfile.write("holepos = ")
    outfile.write(str(holepos))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("J1 = ")
    outfile.write(str(J1[run-runmin]))
    outfile.write("\n")
    outfile.write("J2 = ")
    outfile.write(str(J2[run-runmin]))
    outfile.write("\n")
    outfile.write("J3 = ")
    outfile.write(str(J3[run-runmin]))
    outfile.write("\n")
    outfile.write("J3a = ")
    outfile.write(str(J3a[run-runmin]))
    outfile.write("\n")
    outfile.write("J3b = ")
    outfile.write(str(J3b[run-runmin]))
    outfile.write("\n")
    outfile.write("DXY = ")
    outfile.write(str(DXY[run-runmin]))
    outfile.write("\n")


    outfile.write("ncyc = ")
    outfile.write(str(ncyc))
    outfile.write("\n")
    outfile.write("nequ = ")
    outfile.write(str(nequ))
    outfile.write("\n")
    outfile.write("nbins = ")
    outfile.write(str(nbins))
    outfile.write("\n")
    outfile.write("nconew = ")
    outfile.write(str(nconew))
    outfile.write("\n")
    outfile.write("print = ")
    outfile.write(str(print))
    outfile.write("\n")
    outfile.write("configprint = ")
    outfile.write(str(configprint))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("RANSEED = ")
    outfile.write(str(RANSEED[run-runmin]))
    outfile.write("\n")
    outfile.write("\n")

    outfile.write("EQUENERGY = ")
    outfile.write(str(EQUENERGY))
    outfile.write("\n")
    outfile.write("\n")


    outfile.write("OBC = ")
    outfile.write(str(OBC))
    outfile.write("\n")

    outfile.write("RESETOLDFILES = ")
    outfile.write(str(RESETOLDFILES))
    outfile.write("\n")


for i in range(runmin, runmax+1):
    if mainproject:
        shellfile = open('Jobs/' + mainproject + project + str(i) + '.sh', 'w')
    else:
        shellfile = open('Jobs/' + project + str(i) + '.sh', 'w')
    shellfile.write('#!/bin/bash\n')
    #shellfile.write('valgrind --leak-check=yes ~/Documents/MCHeisenberg/Code/program ' + totalproject + ' ' + str(i) + ' &\n')
    shellfile.write('nice -' + str(NICE) + ' ~/Documents/MCHeisenberg/Code/program ' + totalproject + ' ' + str(i) + ' &\n')
    if mainproject:
        shellfile.write('rm ' + '~/Documents/MCHeisenberg/Jobs/' + mainproject + project + str(i) + '.sh')
    else:
        shellfile.write('rm ' + '~/Documents/MCHeisenberg/Jobs/' + project + str(i) + '.sh')
    shellfile.close()
    if mainproject:
        os.system('chmod +x Jobs/' + mainproject + project + str(i) + '.sh')
    else:
        os.system('chmod +x Jobs/' + project + str(i) + '.sh')
