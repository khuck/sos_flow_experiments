#!/bin/bash  -x
export SOS_CMD_PORT=22500
$* -k 0 &
export SOS_CMD_PORT=22501
$* -k 1 
