#!/bin/bash -x
#PBS -N sos_flow
#PBS -j oe
#PBS -q debug
#PBS -A CSC103
#PBS -l walltime=0:10:00,nodes=257

# set up the environment
. $HOME/src/sos_flow/hosts/ornl/titan/setenv.sh

workdir=$PROJWORK/csc103/khuck/sos_flow_experiments
mkdir -p $workdir
cp $HOME/src/sos_flow/build-titan/bin/sosd $workdir/.
cp $HOME/src/sos_flow/build-titan/bin/showdb $workdir/.
cp $HOME/src/sos_flow_experiments/call_and_response/main $workdir/.
cp $HOME/src/sos_flow_experiments/call_and_response/cleanup.sh $workdir/.
cp $HOME/src/sos_flow_experiments/call_and_response/spawn_aggregators.sh $workdir/.
cd $workdir

hostname=`hostname`
sosbin=$workdir
cwd=`pwd`
num_listeners=256
app_ranks_per_node=2
app_ranks=512
aggregators=2

# export SOS_WORK=${cwd}
if [ -d /dev/shm ] ; then
export SOS_WORK=/dev/shm
else
export SOS_WORK=/tmp
fi

export sos_cmd="${sosbin}/sosd -l ${num_listeners} -a ${aggregators} -w ${SOS_WORK}"
if [ ${num_listeners} == 0 ] ; then
  export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r aggregator"
  export SOS_LISTENER_RANK_OFFSET=0
else
  export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r listener"
  export SOS_LISTENER_RANK_OFFSET=2
fi
export SOS_APP_RANKS_PER_NODE=${app_ranks_per_node}
export SOS_FORK_SHUTDOWN="${sosbin}/sosd_stop"
export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=${cwd}

# clean anything that might be out there
rm -rf ${SOS_WORK}/sosd.* profile* ${SOS_WORK}/start0000*
sleep 2

# launch our aggregator(s)
#aprun -n 1 -N 1 -d 4 ${sos_cmd} -k 0 -r aggregator &
aprun -n 1 -N 1 -d 8 ./spawn_aggregators.sh ${sos_cmd} -r aggregator &
sleep 2

# launch the application
aprun -n ${app_ranks} -N ${app_ranks_per_node} ./main
sleep 2

# get our files
aprun -n ${PBS_NUM_NODES} -N 1 ${cwd}/cleanup.sh
sleep 2

# how did we do?
export SOS_WORK=${cwd}
${sosbin}/showdb
