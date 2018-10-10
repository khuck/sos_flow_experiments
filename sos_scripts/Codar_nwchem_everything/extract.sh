#!/usr/bin/env bash
# Launch the python code that will export the TAU data as an ADIOS BP file.
# Make sure the PYTHONPATH is set to find all the ADIOS, SOS and related python modules.

export SOS_CMD_PORT=22500
export SOS_WORK=`pwd`
export SOS_EVPATH_MEETUP=`pwd`
export SOS_ROOT=/Install/sosflow
export ADIOS_ROOT=/Install/adios-1.13.1
export PYTHONPATH=${SOS_ROOT}/bin:${SOS_ROOT}/lib:${ADIOS_ROOT}/lib/python2.7/site-packages:${PYTHONPATH}
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${SOS_ROOT}/lib:${ADIOS_ROOT}/lib/python2.7/site-packages/adios_mpi:/Install/EVPath/lib

cmd=$*
echo ${cmd}
${cmd} &

# Wait for that to startup
while [ ! -f ${SOS_WORK}/sosd.00000.key ] ; do
    echo "Waiting for ${SOS_WORK}/sosd.00000.key..."
    sleep 1
done

rm -rf ${SOS_WORK}/QA/tests/ethanol/tau-metrics.bp*

echo "Launching ADIOS trace export from SOS..."
python /Codar/nwchem-1/tau_trace_adios.py >& /Codar/nwchem-1/sos_extract.out

# Wait for that to startup
while [ ! -f ${SOS_WORK}/QA/tests/ethanol/tau-metrics.bp ] ; do
    echo "Waiting for ${SOS_WORK}/QA/tests/ethanol/tau-metrics.bp..."
    sleep 1
done
sleep 1

. /etc/profile

echo "Launching Anomaly detection..."
cd /Chimbuko/PerformanceAnalysis/drivers
export PYTHONPATH=/Chimbuko/PerformanceAnalysis/lib/codar/chimbuko/perf_anom:/Install/adios-1.13.1/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
python3 chimbuko.py chimbuko_tau-nwchem.cfg >& /Codar/nwchem-1/anomaly.out

# Once the script has finished, exit the daemon
/Install/sosflow/bin/sosd_stop
