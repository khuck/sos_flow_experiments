#!/bin/bash 
cwd=`pwd`
killall -9 sosd main
rm -rf ${SOS_WORK}/sosd.* profile* ${SOS_WORK}/start0000*
