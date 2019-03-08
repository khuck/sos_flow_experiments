#!/bin/bash -x
set -e

SOS=$HOME/src/sos_flow
rm -f $SOS/sosd.*.db
$SOS/install/bin/sosd -r listener -l 1 -a 0 -k 0 -w $SOS &
cd $HOME/src/tau2/examples/mm
./go.sh