#!/bin/bash -l
#PBS -N test 
#PBS -q debug 
#PBS -A CSC103 
#PBS -l nodes=37,walltime=00:05:00 
#PBS -j oe
#PBS -o both.out

cd ${PBS_O_WORKDIR}
pwd
date

echo "------------------------------------------------------"
echo -n 'Job is running on node '; cat $PBS_NODEFILE
echo "------------------------------------------------------"
echo "PBS: qsub is running on ${PBS_O_HOST}"
echo "PBS: originating queue is ${PBS_O_QUEUE}"
echo "PBS: executing queue is ${PBS_QUEUE}"
echo "PBS: working directory is ${PBS_O_WORKDIR}"
echo "PBS: execution mode is ${PBS_ENVIRONMENT}"
echo "PBS: job identifier is ${PBS_JOBID}"
echo "PBS: job name is ${PBS_JOBNAME}"
echo "PBS: node file is ${PBS_NODEFILE}"
echo "PBS: current home directory is ${PBS_O_HOME}"
echo "PBS: PATH = ${PBS_O_PATH}"
echo "------------------------------------------------------"

export xmainOut="xmain.out" 
export readOut="read2.out"

export SOS_FORK_COMMAND="${PBS_O_WORKDIR}/sosd -l 36 -a 1 -w ${PBS_O_WORKDIR} -k @LISTENER_RANK@ -r listener"
export SOS_FORK_SHUTDOWN="${PBS_O_WORKDIR}/sosd_stop"
export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=${PBS_O_WORKDIR}
export TAU_SOS=1
export TAU_VERBOSE=1

# launch the aggregator
aprun -N 1 -n 1 ./sosd -l 36 -a 1 -w ${PBS_O_WORKDIR} -k 0 -r aggregator &

sleep 5 

# first 32 nodes doing xmain
export SOS_RANKS_PER_NODE=16
export SOS_LISTENER_RANK_OFFSET=1
export PROFILEDIR=profiles_xmain
aprun -N 16 -n 512 ./xmain         > ${xmainOut} 2>&1 &

# last 4 nodes doing reader
export SOS_RANKS_PER_NODE=1
export SOS_LISTENER_RANK_OFFSET=33
export PROFILEDIR=profiles_reader2
aprun -N 1  -n 4   ./reader2 5 x z > ${readOut}  2>&1

# export today=`date | awk '{print $2$3$6}' `
# export logDir=${today}-32m-4r
# mkdir ${logDir}

# rm -f ${logDir}/*

# mv hist* ${logDir}
# mv ${xmainOut} ${logDir}
# mv ${readOut} ${logDir}
# mv both.out ${logDir}
# cp Impact*in ${logDir}

sleep 5
killall -9 aprun

date
