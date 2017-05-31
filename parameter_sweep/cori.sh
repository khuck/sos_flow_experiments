#!/bin/bash
#SBATCH --nodes=2 
#SBATCH -t 5:00:00
#SBATCH -A m1881
#SBATCH -C haswell
#SBATCH --gres=craynetwork:2
#SBATCH -p regular
#######   #SBATCH --qos=interactive

cd ${SCRATCH}/mona/sos_flow_experiments/parameter_sweep
source ${SCRATCH}/mona/sos_flow/hosts/nersc/cori/setenv.sh 

hostname=`hostname`
sosbin=${SCRATCH}/mona/sos_flow/build-cori-icc/bin
cwd=`pwd`
num_listeners=1
app_nodes=1
app_ranks_per_node=32
app_ranks=32

# mytmp=/tmp
mytmp=/dev/shm
#mytmp=${cwd}
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

  # launch the listeners
  if [ ${num_listeners} != 0 ] ; then
    cmd="srun -u -n ${num_listeners} -N ${num_listeners} -c 2 --hint=multithread --nodelist=${othernodes} --gres=craynetwork:1 --mem=25600 ${sos_cmd} -k 1 -r listener"
    echo $cmd
    $cmd &
    sleep 5
  fi
}

stop_servers() {
  # stop the listeners (cascades to aggregator)
  cmd="srun -u -n ${app_nodes} -N ${app_nodes} --nodelist=${othernodes} --gres=craynetwork:1 --mem=25600 ${sosbin}/sosd_stop"
  echo $cmd
  $cmd &
  sleep 5

  # get our files
  srun -u -n 1 -N 1 --gres=craynetwork:1 --mem=25600  --nodelist=${rootnode} ${cwd}/../call_and_response/cleanup.sh
  if [ ${num_listeners} != 0 ] ; then
    srun -u -n ${app_nodes} -N ${app_nodes} --nodelist=${othernodes} --gres=craynetwork:1 --mem=25600 ${cwd}/../call_and_response/cleanup.sh
  fi
  sleep 5

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
cmd_base="--nodelist=${othernodes} ${sosbin}/demo_app"
#ms_per_minute=60000000
ms_per_minute=6000000
#for c in {1..32} ; do
for c in 1 4 8 12 16 20 24 28 32 ; do
#for c in {4..4} ; do
  for p in 1 64 128 192 256 320 384 448 512 ; do
  #for p in 64 ; do
    for d in 1000000 100000 10000 ; do
    #for d in 100000 ; do
      #fresh_start
      launch_servers
      m=$(expr $(expr ${ms_per_minute} / ${d}) \* ${p})
      cmd="srun -u -n ${c} -N ${app_nodes} -c 1 --gres=craynetwork:1 --mem=25600 ${cmd_base} -i 1 -d ${d} -p ${p} -m ${m}"
      echo $cmd
      #${cmd} >& /dev/null
      ${cmd}
      stop_servers ${c} ${p} ${d} ${m}
    done
  done
done
exit
