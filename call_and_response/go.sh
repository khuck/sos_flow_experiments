#!/bin/bash

sosbin=/home/khuck/src/sos_flow/build-linux/bin
cwd=`pwd`
export sos_cmd="${sosbin}/sosd -l 0 -a 1 -w ${cwd}"
export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r aggregator"
export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=${cwd}

killall -9 sosd
rm -rf sosd.* profile* start0000*

# launch the aggregator
# mpirun -np 1 -ppn 1 -hosts n002 ${sos_cmd} -k 0 -r aggregator &
#${sos_cmd} -k 0 -r aggregator &

#sleep 8

# first 32 nodes doing xmain
export SOS_APP_RANKS_PER_NODE=1
export SOS_LISTENER_RANK_OFFSET=0

#mpirun -np 2 ./main
gdb --args ./main
