#!/bin/bash -l
#PBS -N test 
#PBS -q batch 
#PBS -A CSC103 
#PBS -l nodes=33,walltime=00:30:00 
#PBS -j oe
#PBS -o both.out

# echo all commands - for debugging purposes
set -x

###
# -------------------- set the environment -------------------- #
###

# change to the directory where we were launched
cd $PBS_O_WORKDIR

# get the current working directory
export cwd=`pwd`
echo ${cwd}
# timestamp
date

# load modules
module swap PrgEnv-pgi PrgEnv-gnu
module load cudatoolkit
module load papi
module load cmake
module load flexpath/1.12
module load adios/1.12.0
module load python/2.7.9

# set some output file names
export xmainOut="xmain.out" 
export readOut="read2.out"

# SOS_FORK_COMMAND is used when TAU can't connect an SOS client to an SOS listener.
# The application MPI ranks will self-organize, and only 1 rank from each node will
# attempt to launch the SOS listener on that node, using this command.
# @LISTENER_RANK@ will get subsitited for the node rank within all MPI ranks for
# that application, also using an SOS_LISTENER_RANK_OFFSET value. 
export sos_cmd="${cwd}/sosd -l 64 -a 1 -w ${cwd}"
export SOS_FORK_COMMAND="${sos_cmd} -k @LISTENER_RANK@ -r listener"
# Set the TCP port that the listener will listen to, and the port that clients will
# attempt to connect to.
export SOS_CMD_PORT=22500
# Set the directory where the SOS listeners and aggregators will use to establish
# EVPath links to each other.  This needs to be a shared filesystem path, and it
# needs to be writable, of course.
export SOS_EVPATH_MEETUP=${cwd}
# Tell TAU that it should connect to SOS and send TAU data to SOS when adios_close(),
# adios_advance_step() calls are made, and when the application terminates.
export TAU_SOS=1
# If Verbose TAU output is required for debugging.
# export TAU_VERBOSE=1
# Make sure sosd can find libenet.so
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/sw/xk6/flexpath/1.12/cle5.2_gnu4.9.3/lib

# clean up old files, if necessary
rm -rf sosd.* profile* *.out

###
# -------------------- Launch the SOS aggregator -------------- #
###

# launch the aggregator - ALPS will take the first node of the allocation
# and the aggregator will be "rank" 0 in the SOS processes. Launch
# in the background, so we can continue launching other aprun calls.
aprun -n 1 -N 1 ${sos_cmd} -k 0 -r aggregator &

sleep 5 

###
# --------- Launch the listener app in the pipeline ----------- #
###

# 16 nodes doing reader
# ALPS will take the second node of the allocation
# Tell SOS how many application ranks per node there are
export SOS_APP_RANKS_PER_NODE=16
# Tell SOS what "rank" it's listeners should start with - the
# aggregator was "rank" 0, so this node's listener will be 1
export SOS_LISTENER_RANK_OFFSET=1
# Where should TAU write the profile data?
export PROFILEDIR=profiles_reader2
# Go!  in the background, so we can continue launching aprun calls
aprun -n 64 -N $SOS_APP_RANKS_PER_NODE ./reader2 5 x z > ${readOut}  2>&1 &

# Wait a bit.  Just because.
sleep 3

# last 16 nodes doing xmain
# ALPS will take the third, fourth, fifth, sixth nodes of the allocation
# Tell SOS how many application ranks per node there are
export SOS_APP_RANKS_PER_NODE=16
# Tell SOS what "rank" it's listeners should start with - the
# aggregator was "rank" 0, and the reader node was 1, 
# so this node's listeners will start at 33 and be 33,34,35...
export SOS_LISTENER_RANK_OFFSET=17
# Where should TAU write the profile data?
export PROFILEDIR=profiles_xmain
# Go!  in the foreground, so when xmain exits, the PBS allocation will exit.
aprun -n 256 -N 16 ./xmain > ${xmainOut} 2>&1

# wait for clean exit
sleep 3
wait
# timestamp
date
