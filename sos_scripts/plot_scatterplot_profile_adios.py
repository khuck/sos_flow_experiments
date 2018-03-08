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
#from common import TempFile

# Global attr, representing the "previous" frame
previous_attr = dict()
cached_results_dim = None
aggregators = []
sosPort = 0
SOS = None
firstTime = True
cachedName = ""
config = None
last_time = {}
last_fp_ops = {}
last_rss = {}
adios_mode = "w"

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
        print a,sosPort,sql
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

#####
#
def sosScatterplotGenerator():
    global SOS
    global config
    parseConfigFile()
    SOS = SSOS()

    printf("Initializing SOS: ...\b\b\b")
    SOS.init()
    printf("OK!\n")

    # NOTE: When allocation time is scarce, 'stride' here can be
    #       set so that intermediate cycles can be skipped, which is
    #       especially useful when there are thousands of cycles.
    #
    stride = 1

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
    cycleFieldName = "frame"
    #
    num_rows = 0
    # Get at least one active aggregator
    lookupAggregators()
    # wait for a few frames to show up. Frame 0 and 1 are TAU metadata.
    # Frame 2 represents a true iteration.
    max_cycle, maxtime = waitForServer(SOS, cycleFieldName, max(stride,1), True)
    print "Maximum observed '" + cycleFieldName + "' value: " + str(max_cycle) + " (so far)"
    #
    sqlMaxFrame = "SELECT count(*) FROM tblpubs;"
    results, col_names = queryAllAggregators(sqlMaxFrame)
    # We got results from each aggregator, so sum them up
    rank_maxes = [int(x[0]) for x in results]
    rank_max = sum(rank_maxes)
    print "Maximum observed pub_guids: " + str(rank_max)
    #
    #####

    g = None
    if config["output_adios"]:
        # ADIOS file output
        ad.init_noxml()
        g = ad.declare_group("TAU_metrics", "", ad.FLAG.YES)
        ad.define_var(g, "program_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "comm_rank_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "thread_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "metric_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "timer_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "value_count", "", ad.DATATYPE.unsigned_integer, "", "", "")
        ad.define_var(g, "values", "", ad.DATATYPE.unsigned_integer, "value_count,6", "value_count,6", "0,0")
        print "using ADIOS method:", str(config["adios_method"])
        ad.select_method(g, str(config["adios_method"]), "verbose=3", "")

    #
    #
    # EXAMPLE A: Generate .txt set for ALL simulation cycles:
    print "Generating TXT files..."
    lastX = [0.0]*(rank_max + 1)
    lastY = [0.0]*(rank_max + 1)
    lastZ = [0.0]*(rank_max + 1)
    simCycle = max_cycle
    mintime = 0.0
    # Keep running until there are no more frames to wait for.
    # At runtime, this is a moving target, since max_cycle gets updated.
    while config["aggregators"]["runtime"] or simCycle < config["aggregators"]["maxframe"]:
        print "Processing frame", simCycle
        start = time.time()
        vtkOutputFileName = generateADIOSFile(SOS, cycleFieldName, simCycle, lastX, lastY, lastZ, stride, mintime, maxtime, g)
        # clean up the database for long runs
        #cleanDB(SOS, cycleFieldName, simCycle)
        simCycle = simCycle + stride
        # wait for the next batch of frames
        mintime = maxtime
        max_cycle, maxtime = waitForServer(SOS, cycleFieldName, simCycle, False)
        end = time.time()
        print "loop time:", str(end-start)

    #####
    #
    # Whew!  All done!
    #
    # NOTE: See vtkWriter.py for more details.
    #
    if config["output_adios"]:
        ad.finalize()
    SOS.finalize();
    #
    #####
    print "   ...DONE!"
    print 
    return
#
#end:def sosScatterplotGenerator()

def waitForServer(SOS, cycleFieldName, simCycle, first):
    global cachedName
    global config
    # IF running post-processing, skip this "first" - this just makes sure we have data
    # from all publishers.
    if first:
        sum_expected = int(config["aggregators"]["expected_pubs"])
        '''
        sqlFieldNames = "select distinct(name) from tbldata where name like 'TAU..0..calls..%' limit 1;"
        results, col_names = queryAllAggregators(sqlFieldNames)
        cachedName = results[0][0]
        '''
        # cachedName = "TAU..0..calls..MPI"
        cachedName = str(config["events"]["timers"][0])
        # sqlFieldNames = "select count(distinct pub_guid) from viewCombined where value_name = '" + cachedName + "' and " + cycleFieldName + " = " + str(simCycle+1) + ";"
        print "Waiting for publishers", str(simCycle)
        sqlFieldNames = "select count(guid) from tblpubs;"
        waiting = True
        while waiting:
            time.sleep(2.0)
            # query all aggregators
            results, col_names = queryAllAggregators(sqlFieldNames)
            arrived = [int(x[0]) for x in results]
            print arrived, "Publishers have arrived"
            if sum(arrived) >= sum_expected:
                # Wait just a bit more, just in case.
                time.sleep(1.0)
                break

    # how many pubs are there?
    """
    sqlFieldNames = "select count(*) from tblpubs;"
    results, col_names = queryAllAggregators(sqlFieldNames)
    expected = [int(x[0]) for x in results]
    #sum_expected = sum(expected)
    """
    sum_expected = int(config["aggregators"]["expected_pubs"])

    if first or not first:
        # how many pubs have gotten to the expected frame?
        # sqlFieldNames = "select count(pub_guid) from viewCombined where value_name = '" + cachedName + "' and " + cycleFieldName + " = " + str(simCycle) + ";"
        print "Looking for frame", str(simCycle)
        sqlFieldNames = "select count(distinct pub_guid) from viewCombined where " + cycleFieldName + " = " + str(simCycle) + ";"
        maxframe = 0
        waiting = True
        while waiting:
            # query all aggregators
            results, col_names = queryAllAggregators(sqlFieldNames)
            arrived = [int(x[0]) for x in results]
            print arrived
            #if sum(arrived) == sum(expected):
            if sum(arrived) == sum_expected:
                break
            time.sleep(2.0)

    # Everyone has arrived.

    '''
    # What time did the last one finish that frame? Use time_recv, for a little buffer,
    # since we aren't querying *all* variables, just one that we expect.
    sqlFieldNames = "select max(foo) from (select pub_guid, coalesce(max(time_pack),0) as foo from viewCombined where value_name = '" + cachedName + "' and " + cycleFieldName + " = " + str(simCycle) + " group by pub_guid);"
    results, col_names = queryAllAggregators(sqlFieldNames)
    timestamp = [float(x[0]) for x in results]
    maxtime = max(timestamp)
    print "frame",simCycle,"ended at",maxtime
    '''
    maxtime = 0

    return simCycle, maxtime

def cleanDB(SOS, cycleFieldName, simCycle):
    start = time.time()
    sqlstr = "delete from tblvals where " + cycleFieldName + " < " + str(simCycle) + ";"
    results, col_names = queryAllAggregators(sqlstr)
    end = time.time()
    print (end-start), "seconds for cleanup query"

#####
#
def generateADIOSFile(SOS, cycleFieldName, simCycle, lastX, lastY, lastZ, stride, mintime, maxtime, adios_group):
    global previous_attr
    global cached_results_dim
    global cachedName
    global config
    global last_time
    global last_fp_ops
    global last_rss
    global adios_mode

    # Get the frame-specific data...

    # because we are querying different servers, the data has the potential to
    # arrive out-of-order (prog_name, comm_rank). So after each query, sort
    # the results into a dictionary of dictionaries.
    prog_names = {}
    comm_ranks = {}
    value_names = {}
    threads = {}
    metrics = {}
    groups = {}
    timers = {}

    # do the memory first - HWM

    start = time.time()
    sqlValsToColByRank = "select value_name, coalesce(value,0.0), prog_name, comm_rank from viewCombined where value_name like 'TAU_TIMER:%' and frame = " + str(simCycle) + " order by prog_name, comm_rank, value_name;"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    end = time.time()
    print (end-start), "seconds for frame query"

    values_array = np.zeros(shape=(len(results),6), dtype=np.uint32)

    index = 0
    for r in results:
        value_name = str(r[0])
        value = float(r[1])
        prog_name = str(r[2])
        comm_rank = int(r[3])
        if prog_name not in prog_names:
            attr_name = "program_name " + str(len(prog_names))
            prog_names[prog_name] = len(prog_names)
            ad.define_attribute_byvalue(adios_group, attr_name, "", prog_name)
        # may not be necessary...
        if comm_rank not in comm_ranks:
            comm_ranks[comm_rank] = len(comm_ranks)
        # tease apart the timer name.
        tokens = value_name.split(":", 4)
        thread = tokens[1]
        metric = tokens[2]
        group = tokens[3]
        timer = tokens[4]
        if thread not in threads:
            threads[thread] = len(threads)
        if metric not in metrics:
            attr_name = "metric " + str(len(metrics))
            metrics[metric] = len(metrics)
            ad.define_attribute_byvalue(adios_group, attr_name, "", metric)
        if timer not in timers:
            attr_name = "timer " + str(len(timers))
            timers[timer] = len(timers)
            ad.define_attribute_byvalue(adios_group, attr_name, "", timer)
        values_array[index][0] = int(prog_names[prog_name])
        values_array[index][1] = int(comm_ranks[comm_rank])
        values_array[index][2] = int(threads[thread])
        values_array[index][3] = int(metrics[metric])
        values_array[index][4] = timers[timer]
        values_array[index][5] = int(value)
        index = index + 1

    prog_names_array = np.array(list(prog_names.keys()), dtype=np.chararray)
    comm_rank_array = np.array(list(comm_ranks.keys()), dtype=np.uint32)
    thread_array = np.array(list(threads.keys()), dtype=np.uint32)
    metric_array = np.array(list(metrics.keys()), dtype=np.chararray)
    timer_array = np.array(list(timers.keys()), dtype=np.chararray)

    # now that the data is queried and sorted, write it out to the file

    # initialize the ADIOS data
    if config["output_adios"]:
        # write the adios
        fd = ad.open("TAU_metrics", "tau-metrics.bp", adios_mode)
        ad.write_int(fd, "program_count", len(prog_names))
        ad.write_int(fd, "comm_rank_count", len(comm_ranks))
        ad.write_int(fd, "thread_count", len(threads))
        ad.write_int(fd, "metric_count", len(metrics))
        ad.write_int(fd, "timer_count", len(timers))
        ad.write_int(fd, "value_count", len(results))
        ad.write(fd, "values", values_array)
        ad.close(fd)
        #fd.close()
        # future iterations are appending, not writing
        adios_mode = "a"
    """

    return filename
    """
    return ""
#
#end:def generateADIOSFile(...)

def printf(format, *args):
    sys.stdout.write(format % args)

###############################################################################
###############################################################################

if __name__ == "__main__":
    sosScatterplotGenerator()
    #end:FILE


