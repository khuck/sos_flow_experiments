#!/usr/bin/env python

import   time
import   os
import   numpy  as  np
from     ssos   import SSOS
import copy
import glob
import socket
import json

# Global attr, representing the "previous" frame
aggregators = []
sosPort = 0
SOS = None
firstTime = True
config = None

#
# NOTE: This script is intended to be run standalone, and to serve as an
#       example of how to interact with the SOSflow runtime to extract
#       geometry. You may need to update some hardcoded field names that
#       were used in prior experiments for this to work for you.
#

# Method to get the list of running aggregators.
def lookupAggregators():
    global aggregators
    global sosPort
    global firstTime
    global config
    aggregators = []
    # Where are the aggregators located?
    # Port number
    sosPort = os.environ.get("SOS_CMD_PORT")
    # EVPath meetup file
    sosMeetup = os.environ.get("SOS_EVPATH_MEETUP")
    # Wait for the first aggregator, but only for the first query.
    # If the aggregator has disappeared, it may be exiting and we don't
    # want to wait forever
    num_aggregators = int(config["aggregators"]["count"])
    for x in range(num_aggregators):
        first = sosMeetup + "/sosd.0000" + str(x) + ".key"
        print "looking for", first
        while firstTime and not os.path.exists(first):
            time.sleep(1)
    firstTime = False
    # For all key files in meetup location, get the hostname from the key file.
    filter = sosMeetup + "/*.key"
    for f in sorted(glob.glob(filter)):
        fileHandle = open ( f,"r" )
        lineList = fileHandle.readlines()
        fileHandle.close()
        hname = lineList[len(lineList)-1].strip()
        addr = socket.gethostbyname(hname)
        print "Aggregator found:", hname, addr
        aggregators.append(hname)

# Iterate over the running aggregators, and query them.
def queryAllAggregators(sql):
    global sosPort
    global aggregators
    global SOS
    allResults = None
    # Iterate over the aggregators, running the query against each one.
    # This could be done in parallel...
    for a in aggregators:
        #print a,sosPort,sql
        results, col_names = SOS.query(sql, a, sosPort)
        #print "\n...done."
        if allResults == None:
            allResults = results
        else:
            allResults += results
    return allResults, col_names

def parseConfigFile():
    global config
    with open('sos_config.json') as json_data_file:
        config = json.load(json_data_file)
    # set the environment variables
    for var in config["sosd"]:
        if str(var) not in os.environ:
            print "Setting", str(var), "to", config["sosd"][var]
            os.environ[str(var)] = config["sosd"][var]
    # set some defaults
    if "output_text" not in config:
        config["output_text"] = True
    if "output_adios" not in config:
        config["output_adios"] = True
    if "adios_method" not in config:
        config["adios_method"] = "POSIX"
    if "clean_database_after_frame" not in config:
        config["clean_database_after_frame"] = False

def waitForServer(SOS, frame):
    global config

    # how many pubs are there at this frame?
    sum_expected = int(config["aggregators"]["expected_pubs"])

    # how many pubs have gotten to the expected frame?
    print "Looking for frame", str(frame)
    sqlFieldNames = "select count(distinct guid) from tblpubs where latest_frame > " + str(frame) + ";"
    maxframe = 0
    timeouts = 0
    while frame == 0 or timeouts < config["exit_after_n_timeouts"]:
        # query all aggregators
        results, col_names = queryAllAggregators(sqlFieldNames)
        arrived = [int(x[0]) for x in results]
        print arrived, "pubs have arrived at frame", frame
        if sum(arrived) >= sum_expected:
            # Everyone has arrived.
            return False
        time.sleep(config["server_timeout"])
        timeouts = timeouts + 1
    # Too many timeouts, exit.
    return True

def cleanDB(SOS, frame):
    global config
    if config["clean_database_after_frame"]:
        start = time.time()
        sqlstr = "delete from tblvals where frame < " + str(frame) + ";"
        results, col_names = queryAllAggregators(sqlstr)
        end = time.time()
        print (end-start), "seconds for cleanup query"

def do_something(frame):
    start = time.time()
    sqlValsToColByRank = "select prog_name, comm_rank, value_name, coalesce(value,0.0) from viewCombined where value_name like 'TAU_TIMER:%' and frame = " + str(frame) + " order by prog_name, comm_rank, value_name;"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    end = time.time()
    print (end-start), "seconds for query"
    print col_names
    for r in results:
        print r


###############################################################################
###############################################################################

def my_main():
    global SOS
    global config
    parseConfigFile()
    SOS = SSOS()

    print("Initializing SOS: ...")
    SOS.init()
    print("OK!\n")

    # Get at least one active aggregator
    lookupAggregators()
    
    # wait for a frame to show up.
    next_frame = 0
    # while config["aggregators"]["runtime"] or next_frame < config["aggregators"]["maxframe"]:
    while next_frame < config["aggregators"]["maxframe"]:
        # wait for the next batch of frames
        waitForServer(SOS, next_frame)
        print "Processing frame", next_frame
        start = time.time()
        do_something(next_frame)
        # clean up the database for long runs
        cleanDB(SOS, next_frame)
        next_frame = next_frame + 1
        end = time.time()
        print "loop time:", str(end-start)

    # finalize SOS
    SOS.finalize()

    print "   ...DONE!"
    return

if __name__ == "__main__":
    my_main()
