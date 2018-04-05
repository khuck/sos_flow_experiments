#!/bin/bash -l
#PBS -N test 
#PBS -q debug 
#PBS -A CSC103 
#PBS -l nodes=4,walltime=00:10:00 
#PBS -j oe
#PBS -o both.out

set -x

cd $PBS_O_WORKDIR
source /ccs/proj/csc143/khuck/src/sourceme.sh 

export cwd=`pwd`
pwd
date

module swap PrgEnv-pgi PrgEnv-gnu
module load cudatoolkit
module load papi
module load cmake
module load flexpath/1.12
module load adios/1.12.0
module load python/2.7.9

export xmainOut="xmain.out" 
export readOut="read2.out"

export sos_cmd="${cwd}/sosd -l 3 -a 1 -w ${cwd}"
export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r listener"
export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=${cwd}
export TAU_SOS=1
# export TAU_VERBOSE=1
# Make sure sosd can find libenet.so
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/xk6/flexpath/1.12/cle5.2_gnu4.9.3/lib

rm -rf sosd.* profile*

# launch the aggregator
#mpirun -np 1 -ppn 1 -hosts n002 ${sos_cmd} -k 0 -r aggregator &
aprun -n 1 -N 1 ${sos_cmd} -k 0 -r aggregator &

sleep 5 

# last 1 nodes doing reader
export SOS_APP_RANKS_PER_NODE=16
export SOS_LISTENER_RANK_OFFSET=3
export PROFILEDIR=profiles_reader2
aprun -n 16 -N 16 ./reader2 5 x z > ${readOut}  2>&1 &

# first 32 nodes doing xmain
export SOS_APP_RANKS_PER_NODE=32
export SOS_LISTENER_RANK_OFFSET=1
export PROFILEDIR=profiles_xmain
aprun -n 64 -N 32 ./xmain > ${xmainOut} 2>&1 &

date
