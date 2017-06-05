#!/bin/bash
#SBATCH --nodes=3 
#SBATCH -t 0:05:00
#SBATCH -A m1881
#SBATCH -C haswell
#SBATCH -p debug
#######   #SBATCH --gres=craynetwork:2
#######   #SBATCH --qos=interactive

source ${SCRATCH}/mona/sos_flow/hosts/nersc/cori/setenv.sh 
#set -x

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
#export TAU_VERBOSE=1

# Get the first node in the allocation
nodes=`scontrol show hostname $SLURM_NODELIST | paste -d, -s`
IFS=',' read -r rootnode mainnode readernode <<< "$nodes"
echo $rootnode
echo $mainnode
echo $readernode

rm -rf sosd.* profile*

# launch the aggregator
cmd="srun -u -n 1 -N 1 -c 4 --hint=multithread ${sos_cmd} -k 0 -r aggregator"
echo $cmd
$cmd &
sleep 5

# last 4 nodes doing reader
export SOS_APP_RANKS_PER_NODE=1
export SOS_LISTENER_RANK_OFFSET=2
export PROFILEDIR=profiles_reader2
cmd="srun -u -n 1 -N 1 -c 4 --hint=multithread ./reader2 5 x z"
echo ${cmd}
${cmd} > ${readOut} 2>&1 &

# first 32 nodes doing xmain
export SOS_APP_RANKS_PER_NODE=1
export SOS_LISTENER_RANK_OFFSET=1
export PROFILEDIR=profiles_xmain
cmd="srun -u -n 1 -N 1 -c 4 --hint=multithread ./xmain"
echo ${cmd}
${cmd} > ${xmainOut} 2>&1 &

date
