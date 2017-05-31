#!/bin/bash

source $HOME/src/sos_flow/hosts/falcon/setenv.sh 
hostname=`hostname`
sosbin=$ADIOS_ROOT/bin
cwd=`pwd`
num_listeners=1
app_ranks_per_node=12
app_ranks=12
num_nodes=2

unalias rm

export SOS_WORK=${cwd}
#mytmp=/tmp
#mytmp=/dev/shm
#export SOS_WORK=${mytmp}

# Get the first node in the allocation
rootnode=kid91
othernodes=kid90
echo $rootnode
echo $othernodes
lastcore=11

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

fresh_start() {
  # kill anything that might be out there
  cmd="mpirun -np ${num_nodes} --pernode --host kid90,kid91 killall -9 sosd main"
  echo $cmd
  $cmd
  sleep 1

  # clean anything that might be out there
  export SOS_WORK=${mytmp}
  cmd="/bin/rm -rf ${SOS_WORK}/sosd.* profile* ${SOS_WORK}/start0000*"
  echo $cmd
  $cmd
  #sleep 2
}

launch_servers() {
  # launch our aggregator(s)
  cmd="mpirun -np 1 --host ${rootnode} -x SOS_CMD_PORT -x SOS_EVPATH_MEETUP ${sos_cmd} -k 0 -r aggregator"
  echo $cmd
  $cmd &
  sleep 2

  # launch the listeners
  cmd="mpirun -np ${num_listeners} --host ${othernodes} -x SOS_CMD_PORT -x SOS_EVPATH_MEETUP ${sos_cmd} -k 1 -r listener"
  echo $cmd
  $cmd &
  sleep 2
}

stop_servers() {
  # stop the listeners (cascades to aggregator)
  cmd="mpirun -np ${num_listeners} --host ${othernodes} -x SOS_CMD_PORT -x SOS_EVPATH_MEETUP ${sosbin}/sosd_stop"
  echo $cmd
  $cmd &
  sleep 2

  # get our files
  #mpirun -np 2 --pernode --host kid90,kid91 cp ${SOS_WORK}/sosd.* ${cwd} 
  #sleep 2

  export SOS_WORK=${cwd}

     for db in $SOS_WORK/sosd.*.db ; do
       SQL="SELECT 'DATAOUT', '${db}', '${1}', '${2}', '${3}', '${4}', "
       SQL="$SQL COUNT(guid) AS entry_count,"
       SQL="$SQL MIN(time_recv - time_send) AS min_latency,"
       SQL="$SQL AVG(time_recv - time_send) AS avg_latency," 
       SQL="$SQL MAX(time_recv - time_send) AS max_latency, "
       SQL="$SQL AVG(((time_recv - time_send) - sub.a) * ((time_recv - time_send) - sub.a)) AS variance"
       SQL="$SQL FROM tblVals,"
       SQL="$SQL (select avg(time_recv - time_send) as a from tblvals) as sub;"
     #
       #
       sqlite3 -separator ", " $db "$SQL" >> latencies.txt
       # 
    done
}

echo "DATAOUT, 'database', 'ranks', 'pub size', 'delay', 'total size', 'count', 'min latency', 'mean latency', 'max latency', 'variance latency'" > latencies.txt

# launch the applications
# -i pack/publish iterations before delay
# -d delay (milliseconds)
# -p pub size 
# -m total value 
cmd_base="--host ${othernodes} ${sosbin}/demo_app"
#ms_per_minute=60000000
ms_per_minute=6000000
#for c in {1..28} ; do
for c in {1..12} ; do
  for p in 1 64 128 192 256 320 384 448 512 ; do
  #for p in 1 ; do
    for d in 1000000 100000 10000 ; do
    #for d in 1000000 ; do
      fresh_start
      launch_servers
      m=$(expr $(expr ${ms_per_minute} / ${d}) \* ${p})
      cmd="mpirun -np ${c} -x SOS_CMD_PORT -x SOS_EVPATH_MEETUP ${cmd_base} -i 1 -d ${d} -p ${p} -m ${m}"
      #echo $cmd
      time ${cmd} >& /dev/null
      #${cmd}
      stop_servers ${c} ${p} ${d} ${m}
    done
  done
done
exit
