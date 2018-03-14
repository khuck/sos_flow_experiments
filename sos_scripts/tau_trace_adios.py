#!/usr/bin/env python

#####
#
# NOTE: Keep any useful paths to custom Pythons here, with a note.
#
# VPA17 Alpine stack:
# /usr/workspace/wsa/pavis/third_party/toss3_gcc-4.9.3/python/bin/python
#
#####

import   sys
import   subprocess
import   time
import   os
import   random
import   re
import   numpy  as  np
from     ssos   import SSOS
import copy
import glob
import socket
import json
from mpi4py import MPI
import adios_mpi as ad

# Global attr, representing the "previous" frame
aggregators = []
sosPort = 0
SOS = None
firstTime = True
config = None
prog_names = {}
comm_ranks = {}
value_names = {}
threads = {}
metadata_keys = {}
groups = {}
timers = {}
event_types = {}

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
        #sys.stdout.flush()
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

#####
#
def sosToADIOS():
    global SOS
    global config
    parseConfigFile()
    SOS = SSOS()

    printf("Initializing SOS: ...\b\b\b")
    SOS.init()
    printf("OK!\n")

    #####
    #
    # Get the maximum simulation cycle found in the database.
    #
    # NOTE: The cycleFieldName variable should match what is being used
    #       either by your application or SOSflow. If you are not using
    #       an explicit cycle value, you can use SOSflow's internal
    #       field named "frame" that is updated every time SOS_publish(...)
    #       is called. As long as you are publishing to SOS at the end
    #       of major program steps, this will give you what you want.
    #
    # NOTE: For online queries, if you want to ensure that your most
    #       current projection represents a complete set of values,
    #       and you're investigating a block-synchronous code, you can
    #       grab the current maximum and subtract one.
    #
    #
    num_rows = 0
    # Get at least one active aggregator
    lookupAggregators()

    g = None
    if config["output_adios"]:
        # ADIOS file output
        ad.init_noxml()
        g = ad.declare_group("TAU_metrics", "", ad.FLAG.YES)
        ad.define_var(g, "program_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "comm_rank_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "thread_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "event_type_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "timer_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "event_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "event_timestamps", "", ad.DATATYPE.unsigned_long, "event_count,6", "event_count,6", "0,0")
        ad.define_var(g, "comm_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "comm_timestamps", "", ad.DATATYPE.unsigned_long, "event_count,8", "event_count,8", "0,0")
        print "using ADIOS method:", str(config["adios_method"])
        ad.select_method(g, str(config["adios_method"]), "verbose=3", "")

    # wait for a frame to show up. Frame 0 (and maybe 1) are TAU metadata. 
    # The rest should be just timers.
    next_frame = 0
    # first iteration, we are writing the file. after that, appending.
    adios_mode = "w"

    # Keep running until there are no more frames to wait for.
    # At runtime, this is a moving target, since next_frame gets updated.
    while config["aggregators"]["runtime"] or next_frame < config["aggregators"]["maxframe"]:
        # wait for the next batch of frames
        waitForServer(SOS, next_frame)
        print "Processing frame", next_frame
        start = time.time()
        fd = ad.open("TAU_metrics", "tau-metrics.bp", adios_mode)
        writeMetaData(SOS, next_frame, g, fd)
        writeTimerData(SOS, next_frame, g, fd)
        writeCommData(SOS, next_frame, g, fd)
        ad.close(fd)
        # future iterations are appending, not writing
        adios_mode = "a"
        # clean up the database for long runs
        cleanDB(SOS, next_frame)
        next_frame = next_frame + 1
        end = time.time()
        print "loop time:", str(end-start)

    # finalize adios
    if config["output_adios"]:
        ad.finalize()

    # finalize SOS
    SOS.finalize();

    print "   ...DONE!"
    return
#end:def sosToADIOS()

def waitForServer(SOS, frame):
    global config

    # how many pubs are there at this frame?
    sum_expected = int(config["aggregators"]["expected_pubs"])

    # how many pubs have gotten to the expected frame?
    print "Looking for frame", str(frame)
    sqlFieldNames = "select count(distinct guid) from tblpubs where latest_frame > " + str(frame) + ";"
    maxframe = 0
    waiting = True
    while waiting:
        # query all aggregators
        results, col_names = queryAllAggregators(sqlFieldNames)
        arrived = [int(x[0]) for x in results]
        print arrived, "pubs have arrived at frame", frame
        if sum(arrived) >= sum_expected:
            break
        time.sleep(1.0)

    # Everyone has arrived.
    return

def cleanDB(SOS, frame):
    global config
    if config["clean_database_after_frame"]:
        start = time.time()
        sqlstr = "delete from tblvals where frame < " + str(frame) + ";"
        results, col_names = queryAllAggregators(sqlstr)
        end = time.time()
        print (end-start), "seconds for cleanup query"

def writeMetaData(SOS, frame, adios_group, fd):
    global config
    global prog_names
    global comm_ranks
    global value_names
    global threads
    global metadata_keys

    # Get the frame-specific metadata...
    # (there might not be any after frame 0)

    start = time.time()
    sqlValsToColByRank = "select prog_name, comm_rank, value_name, value from viewCombined where value_name like 'TAU_METADATA:%' and frame = " + str(frame) + " order by prog_name, comm_rank, value_name;"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    end = time.time()
    print (end-start), "seconds for frame query"

    for r in results:
        prog_name = str(r[0])
        comm_rank = str(r[1])
        value_name = str(r[2])
        value = str(r[3])
        if prog_name not in prog_names:
            attr_name = "program_name " + str(len(prog_names))
            prog_names[prog_name] = len(prog_names)
            if config["output_adios"]:
                ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, prog_name, "")
        # may not be necessary...
        if comm_rank not in comm_ranks:
            comm_ranks[comm_rank] = len(comm_ranks)
        # tease apart the metadata name.
        tokens = value_name.split(":", 2)
        thread = tokens[1]
        metadata_key = tokens[2]
        if thread not in threads:
            threads[thread] = len(threads)
        attr_name = "MetaData:" + str(prog_names[prog_name]) + ":" + comm_rank + ":" + thread + ":" + metadata_key
        if config["output_adios"]:
            ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, value, "")

    return
#
#end:def writeMetaData(...)

def writeTimerData(SOS, frame, adios_group, fd):
    global config
    global prog_names
    global comm_ranks
    global value_names
    global threads
    global groups
    global timers
    global event_types

    # Get the frame-specific timer data...

    start = time.time()
    sqlValsToColByRank = "select prog_name, comm_rank, value, value_name from viewCombined where (value_name like 'TAU_EVENT_ENTRY%' or value_name like 'TAU_EVENT_EXIT%') and frame = " + str(frame) + " order by value;"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    end = time.time()
    print (end-start), "seconds for frame query"

    values_array = np.zeros(shape=(len(results),6), dtype=np.uint64)

    index = 0
    for r in results:
        prog_name  = str(r[0])
        comm_rank  = str(r[1])
        value     = long(r[2])
        value_name = str(r[3])
        if prog_name not in prog_names:
            attr_name = "program_name " + str(len(prog_names))
            prog_names[prog_name] = len(prog_names)
            ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, prog_name, "")
        # may not be necessary...
        if comm_rank not in comm_ranks:
            comm_ranks[comm_rank] = len(comm_ranks)
        # tease apart the event name
        tokens = value_name.split(":", 2)
        event_type = tokens[0].replace("TAU_EVENT_","")
        if event_type not in event_types:
            attr_name = "event_type " + str(len(event_types))
            event_types[event_type] = len(event_types)
            ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, event_type, "")
        thread = tokens[1]
        timer = tokens[2]
        if thread not in threads:
            threads[thread] = len(threads)
        if timer not in timers:
            attr_name = "timer " + str(len(timers))
            timers[timer] = len(timers)
            ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, timer, "")
        values_array[index][0] = long(prog_names[prog_name])
        values_array[index][1] = long(comm_ranks[comm_rank])
        values_array[index][2] = long(threads[thread])
        values_array[index][3] = long(event_types[event_type])
        values_array[index][4] = long(timers[timer])
        values_array[index][5] = long(value)
        index = index + 1

    # now that the data is queried and in arrays, write it out to the file

    # initialize the ADIOS data
    if config["output_adios"]:
        # write the adios

        # these get written when the comm data is written

        # ad.write_int(fd, "program_count", len(prog_names))
        # ad.write_int(fd, "comm_rank_count", len(comm_ranks))
        # ad.write_int(fd, "thread_count", len(threads))
        #ad.write_int(fd, "timer_count", len(timers))
        #ad.write_int(fd, "event_type_count", len(event_types))

        ad.write_int(fd, "event_count", len(results))
        ad.write(fd, "event_timestamps", values_array)
    return
#
#end:def writeTimerData(...)

def writeCommData(SOS, frame, adios_group, fd):
    global config
    global prog_names
    global comm_ranks
    global value_names
    global threads
    global groups
    global timers
    global event_types

    # Get the frame-specific timer data...

    start = time.time()
    sqlValsToColByRank = "select prog_name, comm_rank, value, value_name from viewCombined where (value_name like 'TAU_EVENT_SEND%' or value_name like 'TAU_EVENT_RECV%') and frame = " + str(frame) + " order by value;"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    end = time.time()
    print (end-start), "seconds for frame query"

    values_array = np.zeros(shape=(len(results),8), dtype=np.uint64)

    index = 0
    for r in results:
        prog_name  = str(r[0])
        comm_rank  = str(r[1])
        value     = long(r[2])
        value_name = str(r[3])
        if prog_name not in prog_names:
            attr_name = "program_name " + str(len(prog_names))
            prog_names[prog_name] = len(prog_names)
            ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, prog_name, "")
        # may not be necessary...
        if comm_rank not in comm_ranks:
            comm_ranks[comm_rank] = len(comm_ranks)
        # tease apart the event name
        tokens = value_name.split(":", 2)
        event_type = tokens[0].replace("TAU_EVENT_","")
        if event_type not in event_types:
            attr_name = "event_type " + str(len(event_types))
            event_types[event_type] = len(event_types)
            ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, event_type, "")
        tokens = value_name.split(":", 4)
        thread = tokens[1]
        tag = tokens[2]
        partner = tokens[3]
        num_bytes = tokens[4]
        if thread not in threads:
            threads[thread] = len(threads)
        values_array[index][0] = long(prog_names[prog_name])
        values_array[index][1] = long(comm_ranks[comm_rank])
        values_array[index][2] = long(threads[thread])
        values_array[index][3] = long(event_types[event_type])
        values_array[index][4] = long(tag)
        values_array[index][5] = long(partner)
        values_array[index][6] = long(num_bytes)
        values_array[index][7] = long(value)
        index = index + 1

    # now that the data is queried and in arrays, write it out to the file

    # initialize the ADIOS data
    if config["output_adios"]:
        # write the adios
        ad.write_int(fd, "program_count", len(prog_names))
        ad.write_int(fd, "comm_rank_count", len(comm_ranks))
        ad.write_int(fd, "thread_count", len(threads))
        ad.write_int(fd, "timer_count", len(timers))
        ad.write_int(fd, "event_type_count", len(event_types))
        ad.write_int(fd, "comm_count", len(results))
        ad.write(fd, "comm_timestamps", values_array)
    return
#
#end:def writeTimerData(...)

def printf(format, *args):
    sys.stdout.write(format % args)

###############################################################################
###############################################################################

if __name__ == "__main__":
    sosToADIOS()
    #end:FILE


