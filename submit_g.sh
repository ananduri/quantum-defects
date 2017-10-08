#!/bin/sh
#PBS -j oe
#PBS -l nodes=1:ppn=8
#PBS -l walltime=23:59:59
#

#OUTDIR=/scratch/network/ananduri/defects



cd $PBS_O_WORKDIR


PROG=getg.py

module load python

python $PROG $N $T $s
