#!/bin/bash -x
#SBATCH --nodes=54 
#SBATCH --sockets-per-node=2 
#SBATCH --cores-per-socket=14 
#SBATCH --partition=defq 

export I_MPI_PMI_LIBRARY=/cm/shared/apps/slurm/16.05.8/lib64/libpmi.so

cd /home/khuck/src/sos_flow_experiments/call_and_response
module load gcc intel slurm
export LD_LIBRARY_PATH=/cm/local/apps/gcc/6.1.0/lib64:$LD_LIBRARY_PATH

source /home/khuck/src/sos_flow/hosts/linux/setenv.sh 
source /home/khuck/src/sos_flow/hosts/linux/hpc.sh 
hostname=`hostname`
sosbin=/home/khuck/src/sos_flow/build-linux/bin
cwd=`pwd`
num_listeners=0
app_ranks_per_node=1
app_ranks=1

# export SOS_WORK=${cwd}
export SOS_WORK=/tmp

# Get the first node in the allocation
nodes=`scontrol show hostname $SLURM_NODELIST | paste -d, -s`
IFS=',' read -r rootnode othernodes <<< "$nodes"
echo $rootnode
echo $othernodes
lastcore=25-27

export sos_cmd="taskset -c ${lastcore} ${sosbin}/sosd -l ${num_listeners} -a 1 -w ${SOS_WORK}"
if [ ${num_listeners} == 0 ] ; then
  export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r aggregator"
  export SOS_LISTENER_RANK_OFFSET=0
else
  export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r listener"
  export SOS_LISTENER_RANK_OFFSET=1
fi
export SOS_APP_RANKS_PER_NODE=${app_ranks_per_node}
export SOS_FORK_SHUTDOWN="${sosbin}/sosd_stop"
export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=${cwd}

# kill anything that might be out there
srun -n 3 -N 3 killall -9 sosd main
sleep 2

# launch our aggregator(s)
srun -n 1 -N 1 --nodelist=${rootnode} ${sos_cmd} -k 0 -r aggregator &
srun -n 1 -N 1 --nodelist=${rootnode} python ${cwd}/example.py
sleep 2

