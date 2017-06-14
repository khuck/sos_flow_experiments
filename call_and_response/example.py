#!/usr/bin/env python
# file "example.py"
#
#   USAGE:   python ./example.py [loop_count]
#
#   See also:   ./trace_example.sh [loop_count]
#
#   Supported SOS.pack(name, type, value) types:
#           SOS.INT
#           SOS.LONG
#           SOS.DOUBLE
#           SOS.STRING
#

import os
import sys
import time
import numpy
import pprint as pp
from ssos import SSOS

last_trigger = 0

def triggerSOSD(SOS):
    print "Analysis: SEND A MESSAGE!"
    sys.stdout.flush()
    sense_handle = "rebalance"
    payload_data = "1"
    payload_size = len(payload_data)
    SOS.trigger(sense_handle, payload_size, payload_data)

def check_balance(SOS, values, iterations):
    global last_trigger
    arr = numpy.array(values)
    themax = numpy.max(arr)
    themean = numpy.mean(arr)
    thedev = numpy.std(arr)
    print themax, themean, thedev
    sys.stdout.flush()
    if thedev > 0.25 * themean and iterations > last_trigger + 2:
        triggerSOSD(SOS)
        last_trigger = iterations
    return themax

def demonstrateSOS():
    SOS = SSOS()

    print "Analysis: Initializing SOS..."
    SOS.init()
    time.sleep(5.0)

    iterations = 0
    while iterations < 48:
        sql_string = "select v.val, v.guid, max(v.time_pack), count(v.time_pack) from tblvals v inner join tbldata d on v.guid = d.guid and d.name like 'Iteration' group by (v.guid);"
        #print "Sending this query to the SOS daemon: "
        #print "    " + sql_string
        results, col_names = SOS.query(sql_string, "localhost", os.environ['SOS_CMD_PORT'])
        values = []
        themax = 5.0
        for i in range(0,len(results)):
            values.append(float(results[i][0]))
            if iterations < int(results[i][3]):
                iterations = int(results[i][3])
        if len(values) > 0:
            themax = check_balance(SOS,values,iterations)
        time.sleep(themax)
   
    SOS.finalize();
    print "Analysis: DONE!"
    print 

if __name__ == "__main__":
    demonstrateSOS()




