#!/bin/bash
#SBATCH --nodes=2 
#SBATCH -t 1:00:00
#SBATCH -A m1881
#SBATCH -C haswell
#SBATCH --gres=craynetwork:2
#SBATCH -p regular
#######   #SBATCH --qos=interactive

cd ${SCRATCH}/call_and_response
source ${SCRATCH}/mona/sos_flow/hosts/nersc/cori/setenv.sh 

hostname=`hostname`
sosbin=${SCRATCH}/mona/sos_flow/build-cori-icc/bin
cwd=`pwd`
num_listeners=1
app_nodes=1
app_ranks_per_node=4
app_ranks=4

# mytmp=/tmp
# mytmp=/dev/shm
mytmp=${cwd}
export SOS_WORK=${mytmp}

# Get the first node in the allocation
nodes=`scontrol show hostname $SLURM_NODELIST | paste -d, -s`
IFS=',' read -r rootnode othernodes <<< "$nodes"
echo $rootnode
echo $othernodes

export sos_cmd="${sosbin}/sosd -l ${num_listeners} -a 1 -w ${SOS_WORK}"
if [ ${num_listeners} == 0 ] ; then
  export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r aggregator"
  export SOS_LISTENER_RANK_OFFSET=0
  othernodes=${rootnode}
else
  export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r listener"
  export SOS_LISTENER_RANK_OFFSET=1
fi
export SOS_APP_RANKS_PER_NODE=${app_ranks_per_node}
export SOS_FORK_SHUTDOWN="${sosbin}/sosd_stop"
export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=${cwd}

fresh_start() {
  # kill anything that might be out there
  cmd="srun -u -n ${SLURM_NNODES} -N ${SLURM_NNODES} killall -9 sosd main"
  echo $cmd
  $cmd
  sleep 2

  # clean anything that might be out there
  export SOS_WORK=${mytmp}
  cmd="srun -u -n ${SLURM_NNODES} -N ${SLURM_NNODES} rm -rf ${SOS_WORK}/sosd.* profile* ${SOS_WORK}/start0000*"
  echo $cmd
  $cmd
  sleep 2
}

launch_servers() {
  # launch our aggregator(s)
  cmd="srun -u -n 1 -N 1 -c 4 --hint=multithread --nodelist=${rootnode} --gres=craynetwork:1 --mem=25600 ${sos_cmd} -k 0 -r aggregator"
  echo $cmd
  $cmd &
  sleep 5
}

fresh_start
launch_servers
# launch the application
cmd="srun -u -n ${app_ranks} -N ${num_listeners} --nodelist=${othernodes} -c 2 --hint=multithread --gres=craynetwork:1 --mem=25600 ./main"
echo $cmd
${cmd}
sleep 2

# get our files
srun -n 1 -N 1 --nodelist=${rootnode} ${cwd}/cleanup.sh
srun -n ${num_listeners} -N ${num_listeners} --nodelist=${othernodes} ${cwd}/cleanup.sh
sleep 2

# how did we do?
export SOS_WORK=${cwd}
${sosbin}/showdb
