#!/bin/bash
#SBATCH --nodes=2 
#SBATCH --sockets-per-node=2 
#SBATCH --cores-per-socket=14 
#SBATCH --partition=defq 

export I_MPI_PMI_LIBRARY=/cm/shared/apps/slurm/16.05.8/lib64/libpmi.so

cd /home/khuck/src/sos_flow_experiments/parameter_sweep
module load gcc intel/17.0.2 slurm
export LD_LIBRARY_PATH=/cm/local/apps/gcc/6.1.0/lib64:$LD_LIBRARY_PATH

source /home/khuck/src/sos_flow/hosts/linux/setenv.sh 
source /home/khuck/src/sos_flow/hosts/linux/hpc.sh 
hostname=`hostname`
sosbin=/home/khuck/src/sos_flow/build-linux/bin
cwd=`pwd`
num_listeners=1
app_ranks_per_node=28
app_ranks=28

# export SOS_WORK=${cwd}
# mytmp=/tmp
mytmp=/dev/shm
export SOS_WORK=${mytmp}

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

fresh_start() {
  # kill anything that might be out there
  cmd="srun -n ${SLURM_NNODES} -N ${SLURM_NNODES} killall -9 sosd main"
  echo $cmd
  $cmd
  sleep 2

  # clean anything that might be out there
  export SOS_WORK=${mytmp}
  cmd="srun -n ${SLURM_NNODES} -N ${SLURM_NNODES} rm -rf ${SOS_WORK}/sosd.* profile* ${SOS_WORK}/start0000*"
  echo $cmd
  $cmd
  sleep 2
}

launch_servers() {
  # launch our aggregator(s)
  cmd="srun -n 1 -N 1 --nodelist=${rootnode} ${sos_cmd} -k 0 -r aggregator"
  echo $cmd
  $cmd &
  sleep 2

  # launch the listeners
  #gperf="--export=LD_PRELOAD=/usr/local/packages/gperftools/2.5/lib/libprofiler.so --export=CPUPROFILE=${cwd}/prof.out"
  export PATH=$PATH:$HOME/src/tau2/x86_64/bin
  gperf="tau_exec -T serial,pthread -ebs"
  cmd="srun -n ${num_listeners} -N ${num_listeners} --nodelist=${othernodes} ${gperf} ${sos_cmd} -k 1 -r listener"
  echo $cmd
  $cmd &
  sleep 2
}

stop_servers() {
  # stop the listeners (cascades to aggregator)
  cmd="srun -n ${num_listeners} -N ${num_listeners} --nodelist=${othernodes} ${sosbin}/sosd_stop"
  echo $cmd
  $cmd &
  sleep 2

  # get our files
  srun -n 1 -N 1 --nodelist=${rootnode} ${cwd}/../call_and_response/cleanup.sh
  srun -n ${num_listeners} -N ${num_listeners} --nodelist=${othernodes} ${cwd}/../call_and_response/cleanup.sh
  sleep 2

  export SOS_WORK=${cwd}

     for db in $SOS_WORK/sosd.*.db ; do
       SQL="SELECT 'DATAOUT', '${db}', '${1}', '${2}', '${3}', '${4}', "
       SQL="$SQL MAX(rowid) AS entry_count,"
       SQL="$SQL MIN(time_recv - time_send) AS min_latency,"
       SQL="$SQL AVG(time_recv - time_send) AS avg_latency," 
       SQL="$SQL MAX(time_recv - time_send) AS max_latency "
       SQL="$SQL FROM tblVals;"
     #
       #
       sqlite3 -separator ", " $db "$SQL" >> latencies.txt
       # 
    done
}

echo "DATAOUT, 'database', 'ranks', 'pub size', 'delay', 'total size', 'min latency', 'mean latency', 'max latency'" > latencies.txt

# launch the applications
# -i pack/publish iterations before delay
# -d delay (milliseconds)
# -p pub size 
# -m total value 
cmd_base="--nodelist=${othernodes} ${sosbin}/demo_app"
#ms_per_minute=60000000
ms_per_minute=6000000
#for c in {1..28} ; do
for c in {9..9} ; do
  #for p in 1 64 1024 16384 ; do
  for p in 16384 ; do
    #for d in 1000000 100000 10000 ; do
    for d in 10000 ; do
      fresh_start
      launch_servers
      m=$(expr $(expr ${ms_per_minute} / ${d}) \* ${p})
      cmd="srun -n ${c} -N ${num_listeners} ${cmd_base} -i 1 -d ${d} -p ${p} -m ${m}"
      echo $cmd
      ${cmd} >& /dev/null
      stop_servers ${c} ${p} ${d} ${m}
    done
  done
done
exit
