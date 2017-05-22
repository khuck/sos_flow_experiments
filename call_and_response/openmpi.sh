#!/bin/bash 

hostname=`hostname`
sosbin=/home/khuck/src/sos_flow/build-linux/bin
cwd=`pwd`
num_listeners=0
app_ranks_per_node=4
app_ranks=4

export sos_cmd="taskset -c 7 ${sosbin}/sosd -l ${num_listeners} -a 1 -w ${cwd}"
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

killall -9 sosd mpirun main python
rm -rf sosd.* profile* start0000*

# launch the aggregator
# ${sos_cmd} -k 0 -r aggregator &
#mpirun -np 1 --host n002 ${sos_cmd} -k 0 -r aggregator &

sleep 2

mpirun -np ${app_ranks} --map-by core ./main &

sleep 2
python ./example.py

# ${sosbin}/showdb
