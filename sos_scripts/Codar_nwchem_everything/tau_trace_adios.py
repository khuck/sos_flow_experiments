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
from collections import Counter
import json
from mpi4py import MPI
import adios_mpi as ad
from operator import itemgetter

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
counters = {}
event_types = {}
column_map = {}

validation = {}

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
        print ("looking for", first)
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
        print ("Aggregator found:", hname, addr)
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

def queryAllAggregatorsCache(pub_filter, val_filter, frame_start, frame_depth):
    global sosPort
    global aggregators
    global SOS
    allResults = None
    # Iterate over the aggregators, running the query against each one.
    # This could be done in parallel...
    for a in aggregators:
        print (a,sosPort,"pub: '",pub_filter,"' val: '",val_filter,"' frame:",frame_start," depth:",frame_depth)
        results, col_names = SOS.cache_grab(pub_filter,val_filter,frame_start,frame_depth,a, sosPort)
        #print "\n...done."
        if allResults == None:
            allResults = results
        else:
            allResults += results
    return allResults, col_names

def queryAllAggregatorsManifest(pub_filter):
    global sosPort
    global aggregators
    global SOS
    allResults = None
    # Iterate over the aggregators, running the query against each one.
    # This could be done in parallel...
    frames = []
    for a in aggregators:
        max_frame, results, col_names = SOS.request_pub_manifest(pub_filter,a,sosPort)
        #print "\n...done."
        if allResults == None:
            allResults = results
        else:
            allResults += results
        frames.append(max_frame)
    return allResults, col_names, min(frames)


def parseConfigFile():
    global config
    with open('sos_config.json') as json_data_file:
        config = json.load(json_data_file)
    # set the environment variables
    for var in config["sosd"]:
        if str(var) not in os.environ:
            print ("Setting", str(var), "to", config["sosd"][var])
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
    if "server_timeout" not in config:
        # Specified in seconds
        config["server_timeout"] = 1
    if "exit_after_n_timeouts" not in config:
        config["exit_after_n_timeouts"] = 100

#####
#
def sosToADIOS():
    global SOS
    global config
    global validation
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
        ad.define_var(g, "timer_event_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "event_timestamps", "", ad.DATATYPE.unsigned_long, "timer_event_count,6", "timer_event_count,6", "0,0")
        ad.define_var(g, "counter_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "counter_event_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "counter_values", "", ad.DATATYPE.unsigned_long, "counter_event_count,6", "counter_event_count,6", "0,0")
        ad.define_var(g, "comm_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "comm_timestamps", "", ad.DATATYPE.unsigned_long, "comm_count,8", "comm_count,8", "0,0")
        print ("using ADIOS method:", str(config["adios_method"]))
        ad.select_method(g, str(config["adios_method"]), "verbose=3", "")

    # wait for a frame to show up. Frame 0 (and maybe 1) are TAU metadata. 
    # The rest should be just timers.
    next_frame = 0
    # first iteration, we are writing the file. after that, appending.
    adios_mode = "w"

    waitForServer(SOS, 0)
    buildColumnMap(SOS)

    # Keep running until there are no more frames to wait for.
    # At runtime, this is a moving target, since next_frame gets updated.
    done = False
    total_count = 0
    while (not done or total_count > 0) and (config["aggregators"]["runtime"] or next_frame < config["aggregators"]["maxframe"]):
        # wait for the next batch of frames
        if not done:
            timeout = waitForServer(SOS, next_frame + 1)
        if timeout:
            done = True
        #if len(column_map) == 0:
        #    buildColumnMap(SOS)
        print ("Processing frame", next_frame)
        start = time.time()
        fd = ad.open("TAU_metrics", "tau-metrics.bp", adios_mode)
        meta_count = writeMetaData(SOS, next_frame, g, fd)
        timer_count = writeTimerData(SOS, next_frame, g, fd)
        total_count = meta_count + timer_count
        ad.close(fd)
        # future iterations are appending, not writing
        adios_mode = "a"
        print ("Processed", total_count, "rows")
        if total_count == 0 and done:
            break
        next_frame = next_frame + 1
        end = time.time()
        print ("loop time:", str(end-start))

    # finalize adios
    if config["output_adios"]:
        ad.finalize()

    # finalize SOS
    SOS.finalize();

    for p in validation:
        for r in validation[p]:
            for t in validation[p][r]:
                if len(validation[p][r][t]) != 0:
                    print ("VALIDATION ERROR!", p, r, t, validation[p][r][t], "was not exited")
    print ("   ...DONE!")
    return
#end:def sosToADIOS()

def buildColumnMap(SOS):
    global config
    global column_map

    # how many pubs are there at this frame?
    sum_expected = int(config["aggregators"]["expected_pubs"])

    pub_filter = ""
    val_filter = ""
    frame_start = -1
    frame_depth = -1

    print ("Getting column indices")
    col_names = []
    while len(col_names) == 0:
        # query all aggregators
        results, col_names = queryAllAggregatorsCache(pub_filter, val_filter, frame_start, frame_depth)
        time.sleep(config["server_timeout"])

    column_map = {}
    for c in col_names:
        column_map[str(c)] = len(column_map)

def waitForServer(SOS, frame):
    global config

    # how many pubs are there at this frame?
    sum_expected = int(config["aggregators"]["expected_pubs"])

    pub_filter = ""
    val_filter = ""
    frame_start = frame
    frame_depth = 1

    # how many pubs have gotten to the expected frame?
    print ("Looking for frame", str(frame))
    maxframe = 0
    timeouts = 0
    while frame <= 1 or timeouts < config["exit_after_n_timeouts"]:
        # query all aggregators
        results, col_names, max_frame = queryAllAggregatorsManifest(pub_filter)
        if len(results) > 0:
            print (col_names)
            frame_column = col_names.index("pub_frame")
            frames = [int(x[frame_column]) for x in results]
            # How many unique pub_guid values do we have?
            count = Counter(frames)
            arrived = 0
            for f in count:
                if f >= frame:
                    arrived = arrived + count[f]
            print (arrived, "pubs have arrived at frame", frame)
            if arrived >= sum_expected:
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
        print ((end-start), "seconds for cleanup query")

def writeMetaData(SOS, frame, adios_group, fd):
    global config
    global prog_names
    global comm_ranks
    global value_names
    global threads
    global metadata_keys
    global column_map

    # Get the frame-specific metadata...
    # (there might not be any after frame 0)

    start = time.time()
    # sqlValsToColByRank = "select prog_name, comm_rank, value_name, value from viewCombined where value_name like 'TAU_METADATA:%' and frame = " + str(frame) + " order by prog_name, comm_rank, value_name;"
    #results, col_names = queryAllAggregators(sqlValsToColByRank)

    pub_filter = ""
    val_filter = "TAU_Metadata"
    frame_start = frame
    frame_depth = 1
    results, col_names = queryAllAggregatorsCache(pub_filter, val_filter, frame_start, frame_depth)

    end = time.time()
    print ((end-start), "seconds for metadata query")
    prog_name_index = column_map["prog_name"]
    comm_rank_index = column_map["comm_rank"]
    value_name_index = column_map["val_name"]
    value_index = column_map["val"]
    frame_index = column_map["frame"]
    total_valid = len(results)

    for r in results:
        prog_name = str(r[prog_name_index])
        comm_rank = str(r[comm_rank_index])
        value_name = str(r[value_name_index])
        value = str(r[value_index])
        this_frame = int(r[frame_index])
        if this_frame != frame:
            total_valid = total_valid - 1
            continue
        if value == "":
            #print "skipping", value_name
            continue
        if prog_name not in prog_names:
            attr_name = "program_name " + str(len(prog_names))
            prog_names[prog_name] = len(prog_names)
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
        #print attr_name,value
        ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, value, "")

    return total_valid
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
    global validation

    # Get the frame-specific timer data...

    start = time.time()
    #sqlValsToColByRank = "select prog_name, comm_rank, value, value_name from viewCombined where (value_name like 'TAU_EVENT_ENTRY%' or value_name like 'TAU_EVENT_EXIT%') and frame = " + str(frame) + " order by value;"
    #results, col_names = queryAllAggregators(sqlValsToColByRank)
    pub_filter = ""
    val_filter = "TAU_EVENT"
    frame_start = frame
    frame_depth = 1
    results, col_names = queryAllAggregatorsCache(pub_filter, val_filter, frame_start, frame_depth)

    end = time.time()
    print ((end-start), "seconds for event query")

    timer_values_array = np.zeros(shape=(len(results),6), dtype=np.uint64)
    counter_values_array = np.zeros(shape=(len(results),6), dtype=np.uint64)
    comm_values_array = np.zeros(shape=(len(results),8), dtype=np.uint64)

    timer_index = 0
    counter_index = 0
    comm_index = 0
    prog_name_index = column_map["prog_name"]
    comm_rank_index = column_map["comm_rank"]
    value_name_index = column_map["val_name"]
    value_index = column_map["val"]
    frame_index = column_map["frame"]
    time_index = column_map["time_pack"]
    results = sorted(results, key=itemgetter(value_index))
    total_valid = len(results)

    for r in results:
        prog_name  = str(r[prog_name_index])
        comm_rank  = str(r[comm_rank_index])
        value     = int(r[value_index])
        value_name = str(r[value_name_index])
        row_frame = str(r[frame_index])
        #print row_frame, prog_name, comm_rank, value_name
        if int(row_frame) != frame:
            total_valid = total_valid - 1
            continue
        if prog_name not in validation:
            validation[prog_name] = {}
        if comm_rank not in validation[prog_name]:
            validation[prog_name][comm_rank] = {}
        if "TAU_EVENT_ENTRY" in value_name or "TAU_EVENT_EXIT" in value_name:
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
            thread = int(tokens[1])
            if thread not in validation[prog_name][comm_rank]:
                validation[prog_name][comm_rank][thread] = []
            timer = tokens[2]
            if thread not in threads:
                threads[thread] = len(threads)
            if timer not in timers:
                attr_name = "timer " + str(len(timers))
                timers[timer] = len(timers)
                ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, timer, "")
            """
            if "MPI_Send" in value_name:
                print (value_name, thread)
            if "MPI_Recv" in value_name:
                print (value_name, thread)
            """
            timer_values_array[timer_index][0] = int(prog_names[prog_name])
            timer_values_array[timer_index][1] = int(comm_ranks[comm_rank])
            timer_values_array[timer_index][2] = int(thread)
            timer_values_array[timer_index][3] = int(event_types[event_type])
            timer_values_array[timer_index][4] = int(timers[timer])
            timer_values_array[timer_index][5] = int(value)
            timer_index = timer_index + 1
            if "TAU_EVENT_ENTRY" in value_name:
                validation[prog_name][comm_rank][thread].append(timer)
            else:
                if len(validation[prog_name][comm_rank][thread]) == 0:
                    print ("VALIDATION ERROR! empty stack", prog_name, comm_rank, thread, timer)
                    #sys.exit()
                else:
                    current_timer = validation[prog_name][comm_rank][thread].pop()
                    if current_timer != timer:
                        print ("VALIDATION ERROR!", value, prog_names[prog_name], comm_rank, thread, timers[timer], "!= current: ", timers[current_timer])
        elif "TAU_EVENT_COUNTER" in value_name:
            # convert the timestamp from seconds to usec
            timestamp = float(r[time_index]) * 1000000
            if prog_name not in prog_names:
                attr_name = "program_name " + str(len(prog_names))
                prog_names[prog_name] = len(prog_names)
                ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, prog_name, "")
            # may not be necessary...
            if comm_rank not in comm_ranks:
                comm_ranks[comm_rank] = len(comm_ranks)
            # tease apart the event name
            tokens = value_name.split(":", 2)
            thread = tokens[1]
            counter = tokens[2]
            if thread not in threads:
                threads[thread] = len(threads)
            if counter not in counters:
                attr_name = "counter " + str(len(counters))
                counters[counter] = len(counters)
                ad.define_attribute(adios_group, attr_name, "", ad.DATATYPE.string, counter, "")
            counter_values_array[counter_index][0] = int(prog_names[prog_name])
            counter_values_array[counter_index][1] = int(comm_ranks[comm_rank])
            counter_values_array[counter_index][2] = int(thread)
            counter_values_array[counter_index][3] = int(counters[counter])
            counter_values_array[counter_index][4] = int(value)
            counter_values_array[counter_index][5] = int(timestamp)
            counter_index = counter_index + 1
        elif "TAU_EVENT_SEND" in value_name or "TAU_EVENT_RECV" in value_name:
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
            comm_values_array[comm_index][0] = int(prog_names[prog_name])
            comm_values_array[comm_index][1] = int(comm_ranks[comm_rank])
            comm_values_array[comm_index][2] = int(thread)
            comm_values_array[comm_index][3] = int(event_types[event_type])
            comm_values_array[comm_index][4] = int(tag)
            comm_values_array[comm_index][5] = int(partner)
            comm_values_array[comm_index][6] = int(num_bytes)
            comm_values_array[comm_index][7] = int(value)
            comm_index = comm_index + 1
        else:
            print ("ERROR! unknown event:", prog_name, comm_rank, value_name)
    # now that the data is queried and in arrays, write it out to the file

    # initialize the ADIOS data
    if config["output_adios"] and (timer_index > 0 or counter_index > 0 or comm_index > 0):
    #if config["output_adios"]:
        # write the adios
        ad.write_int(fd, "program_count", len(prog_names))
        ad.write_int(fd, "comm_rank_count", len(comm_ranks))
        ad.write_int(fd, "thread_count", len(threads))
        ad.write_int(fd, "timer_count", len(timers))
        ad.write_int(fd, "event_type_count", len(event_types))
        ad.write_int(fd, "timer_event_count", timer_index)
        ad.write_int(fd, "counter_count", len(counters))
        ad.write_int(fd, "counter_event_count", counter_index)
        ad.write_int(fd, "comm_count", comm_index)
        np.resize(timer_values_array, (timer_index,6))
        np.resize(counter_values_array, (counter_index,6))
        np.resize(comm_values_array, (comm_index,8))
        ad.write(fd, "event_timestamps", timer_values_array)
        ad.write(fd, "counter_values", counter_values_array)
        ad.write(fd, "comm_timestamps", comm_values_array)
    return total_valid
#
#end:def writeTimerData(...)

def printf(format, *args):
    sys.stdout.write(format % args)

###############################################################################
###############################################################################

if __name__ == "__main__":
    sosToADIOS()

    #end:FILE


