#!/bin/bash -l
#PBS -N test 
#PBS -q debug 
#PBS -A CSC103 
#PBS -l nodes=3,walltime=00:10:00 
#PBS -j oe
#PBS -o both.out

set -x

export cwd=`pwd`
pwd
date

export xmainOut="xmain.out" 
export readOut="read2.out"

export sos_cmd="${cwd}/sosd -l 2 -a 1 -w ${cwd}"
export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r listener"
export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=${cwd}
export TAU_SOS=1
export TAU_VERBOSE=1
export TAU_CONFIG=-gnu-mpi-pthread
export TAU_ARCH=x86_64
export TAU_ROOT=/home/khuck/src/tau2
export TAU_MAKEFILE=${TAU_ROOT}/${TAU_ARCH}/lib/Makefile.tau${TAU_CONFIG}
export PATH=${TAU_ROOT}/${TAU_ARCH}/bin:${PATH}
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${HOME}/src/chaos/linux-gcc/lib

rm -rf sosd.* profile*

# launch the aggregator
#mpirun -np 1 -ppn 1 -hosts n002 ${sos_cmd} -k 0 -r aggregator &
${sos_cmd} -k 0 -r aggregator &

sleep 5 

tauexec="tau_exec -v -T pdt,mpi,pthread"
#tauexec=""

# last 4 nodes doing reader
export SOS_APP_RANKS_PER_NODE=1
export SOS_LISTENER_RANK_OFFSET=2
export PROFILEDIR=profiles_reader2
mpirun -np 1 -ppn 1 -hosts n004 ${tauexec} ./reader2 5 x z > ${readOut}  2>&1 &

# first 32 nodes doing xmain
export SOS_APP_RANKS_PER_NODE=1
export SOS_LISTENER_RANK_OFFSET=1
export PROFILEDIR=profiles_xmain
mpirun -np 1 -ppn 1 -hosts n003 ${tauexec} ./xmain         > ${xmainOut} 2>&1 &
#mpirun -np 1 -ppn 1 -hosts n003 ${tauexec} gdb --args ./xmain 
#aprun -N 1 -n 1 ./xmain &

# export today=`date | awk '{print $2$3$6}' `
# export logDir=${today}-32m-4r
# mkdir ${logDir}

# rm -f ${logDir}/*

# mv hist* ${logDir}
# mv ${xmainOut} ${logDir}
# mv ${readOut} ${logDir}
# mv both.out ${logDir}
# cp Impact*in ${logDir}

# sleep 5
# killall -9 aprun

date
