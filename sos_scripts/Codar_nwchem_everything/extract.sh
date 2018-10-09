#!/usr/bin/env bash
# Launch the python code that will export the TAU data as an ADIOS BP file.
# Make sure the PYTHONPATH is set to find all the ADIOS, SOS and related python modules.

cmd=$*
echo ${cmd}
${cmd} &

export SOS_CMD_PORT=22500
export SOS_WORK=`pwd`
export SOS_EVPATH_MEETUP=`pwd`

# Wait for that to startup
while [ ! -f ${SOS_WORK}/sosd.00000.key ] ; do
    echo "Waiting for ${SOS_WORK}/sosd.00000.key..."
    sleep 1
done

export SOS_ROOT=/Install/sosflow
export ADIOS_ROOT=/Install/adios-1.13.1
export PYTHONPATH=${SOS_ROOT}/bin:${SOS_ROOT}/lib:${ADIOS_ROOT}/lib/python2.7/site-packages:${PYTHONPATH}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SOS_ROOT}/lib:${ADIOS_ROOT}/lib/python2.7/site-packages/adios_mpi:/Install/EVPath/lib

echo "Launching ADIOS trace export from SOS..."
#python ./tau_trace_adios.py
python ./tau_trace_adios_cache.py
