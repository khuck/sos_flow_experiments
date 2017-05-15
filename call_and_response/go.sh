#!/bin/bash -x
#SBATCH --nodes=3 
#SBATCH --sockets-per-node=2 
#SBATCH --cores-per-socket=14 
#SBATCH --partition=defq 

export I_MPI_PMI_LIBRARY=/cm/shared/apps/slurm/16.05.8/lib64/libpmi.so

cd /home/khuck/src/sos_flow_experiments/call_and_response
module load gcc intel slurm
export LD_LIBRARY_PATH=/cm/local/apps/gcc/6.1.0/lib64:$LD_LIBRARY_PATH

hostname=`hostname`
sosbin=/home/khuck/src/sos_flow/build-linux/bin
cwd=`pwd`
num_listeners=2
app_ranks_per_node=1
app_ranks=2

export sos_cmd="${sosbin}/sosd -l ${num_listeners} -a 1 -w ${cwd}"
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

killall -9 sosd
rm -rf sosd.* profile* start0000*

# launch the aggregator
# ${sos_cmd} -k 0 -r aggregator &
#mpirun -np 1 --host n002 ${sos_cmd} -k 0 -r aggregator &
srun -n 1 -N 1 -c 28 --nodelist=${SLURMD_NODENAME} ${sos_cmd} -k 0 -r aggregator &
#${sos_cmd} -k 0 -r aggregator &

sleep 2

#mpirun -np ${app_ranks} -ppn ${app_ranks_per_node} -hosts n003,n004 ./main
#srun -n ${app_ranks} -c ${cpus_per_task} -nodelist=n003,n004 ./main
srun -n ${app_ranks} -N ${app_ranks_per_node} --exclude=${SLURMD_NODENAME} ./main
#gdb --args ./main
#./main

#mpirun -np 3 -ppn 1 killall -9 sosd
#mpirun -np 3 -ppn 1 ${sosbin}/sosd_stop
${sosbin}/showdb
