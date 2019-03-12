#!/bin/bash -x
set -e

SOS=$HOME/src/sos_flow
rm -f $SOS/sosd.*.db

export SOS_CMD_PORT=22500
export SOS_EVPATH_MEETUP=$SOS
export SOS_DB_DISABLED=TRUE
export SOS_PUB_CACHE_DEPTH=1000
$SOS/install/bin/sosd -r listener -l 1 -a 0 -k 0 -w $SOS &

echo "Waiting on $SOS/sosd.00000.key..."
while [ ! -f $SOS/sosd.00000.key ]
do
  sleep 2
done

cd $HOME/src/tau2/examples/mm
./go.sh