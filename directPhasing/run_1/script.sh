#!/bin/bash
#
#PBS -N 1no4
#PBS -o xxx.log
#PBS -e xxx.err
#PBS -l mem=120GB
#PBS -l walltime=99:00:00

export FI_TCP_IFACE=em2 #make em2 as the only connection between the other nodes

# Calculate the number of processors allocated to this run.
NPROCS=`wc -l < $PBS_NODEFILE`

# Calculate the number of nodes allocated.
NNODES=`uniq $PBS_NODEFILE | wc -l`

### Display the job context
echo Running on host `hostname`
echo Time is `date`


cd $PBS_O_WORKDIR
mpiicpc *.cpp -mcmodel=large -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -O -no-multibyte-chars;
mpirun -np 100 -f $PBS_NODEFILE ./a.out
