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
#import   numpy  as  np
from     ssos   import SSOS
import copy
import glob
import socket
import json
from threading import Thread

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

def threadedQuery(a,sosPort,sql,index,tmp_results):
    print a,sosPort,sql,
    results, col_names = SOS.query(sql, a, sosPort)
    print "thread", str(index), "...done."
    tmp_results[index] += results

# Iterate over the running aggregators, and query them.
def queryAllAggregators(sql):
    global sosPort
    global aggregators
    global SOS
    allResults = None
    # set up our temporary results
    num_aggregators = len(aggregators)
    tmp_results = []
    for x in range(num_aggregators):
        tmp_results.append([])

    # Iterate over the aggregators, running the query against each one.
    # This could be done in parallel...
    index = 0
    threads = []
    for a in aggregators:
        process = Thread(target=threadedQuery, args=[a,sosPort,sql,index,tmp_results])
        process.start()
        threads.append(process)
        index = index + 1

    index = 0
    allResults = []
    for process in threads:
        process.join()
        allResults += tmp_results[index]
        index = index + 1

    col_names = [] # not used
    return allResults, col_names

def parseConfigFile():
    global config
    with open('sos_config.json') as json_data_file:
        config = json.load(json_data_file)
    # set the environment variables
    for var in config["sosd"]:
        os.environ[str(var)] = config["sosd"][var]

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
        printf("Processing frame %d ", simCycle)
        start = time.time()
        vtkOutputFileName = generateVTKFile(SOS, cycleFieldName, simCycle, lastX, lastY, lastZ, stride, mintime, maxtime)
        simCycle = simCycle + stride
        #if simCycle == 100:
        #    break;
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
        '''
        sqlFieldNames = "select distinct(name) from tbldata where name like 'TAU::0::calls::%' limit 1;"
        results, col_names = queryAllAggregators(sqlFieldNames)
        cachedName = results[0][0]
        '''
        # cachedName = "TAU::0::calls::MPI"
        cachedName = str(config["events"]["timers"][0])
        # sqlFieldNames = "select count(distinct pub_guid) from viewCombined where value_name = '" + cachedName + "' and " + cycleFieldName + " = " + str(simCycle+1) + ";"
        print "Looking for frame", str(simCycle)
        sqlFieldNames = "select count(distinct pub_guid) from viewCombined where " + cycleFieldName + " = " + str(simCycle) + ";"
        waiting = True
        while waiting:
            # query all aggregators
            results, col_names = queryAllAggregators(sqlFieldNames)
            arrived = [int(x[0]) for x in results]
            time.sleep(1.0)
            if sum(arrived) > 0:
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
            time.sleep(1.0)

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

#####
#
def generateVTKFile(SOS, cycleFieldName, simCycle, lastX, lastY, lastZ, stride, mintime, maxtime):
    global previous_attr
    global cached_results_dim
    global cachedName
    global config
    global last_time
    global last_fp_ops

    # Get the frame-specific data...

    # because we are querying different servers, the data has the potential to
    # arrive out-of-order (prog_name, comm_rank). So after each query, sort
    # the results into a dictionary of dictionaries.
    prog_names = {}

    # do the memory first

    start = time.time()
    sqlValsToColByRank = "select coalesce(value,0.0), prog_name, comm_rank from viewCombined where value_name = '" + str(config["events"]["counters"][0]) + "' and frame = " + str(simCycle) + ";"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    print "query done"
    end = time.time()
    print (end-start), "seconds for frame query"

    for r in results:
        value = float(r[0])
        prog_name = str(r[1])
        comm_rank = int(r[2])
        if prog_name not in prog_names:
            prog_names[prog_name] = {}
        if comm_rank not in prog_names[prog_name]:
            prog_names[prog_name][comm_rank] = []
        prog_names[prog_name][comm_rank].append(value)

    # do the timer PAPI metric next

    start = time.time()
    # Make sure we sort by prog_name, comm_rank!
    sqlValsToColByRank = "select coalesce(value,0.0), prog_name, comm_rank from viewCombined where value_name = '" + str(config["events"]["timers"][0]) + "' and frame = " + str(simCycle) + ";"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    print "query done"
    end = time.time()
    print (end-start), "seconds for frame query"

    for r in results:
        value = float(r[0])
        prog_name = str(r[1])
        comm_rank = int(r[2])
        if prog_name not in prog_names:
            prog_names[prog_name] = {}
        if comm_rank not in prog_names[prog_name]:
            prog_names[prog_name][comm_rank] = []
        prog_names[prog_name][comm_rank].append(value)

    # do the timer time metric next

    start = time.time()
    sqlValsToColByRank = "select coalesce(value,0.0), prog_name, comm_rank from viewCombined where value_name = '" + str(config["events"]["timers"][1]) + "' and frame = " + str(simCycle) + ";"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    print "query done"
    end = time.time()
    print (end-start), "seconds for frame query"

    for r in results:
        value = float(r[0])
        prog_name = str(r[1])
        comm_rank = int(r[2])
        if prog_name not in prog_names:
            prog_names[prog_name] = {}
        if comm_rank not in prog_names[prog_name]:
            prog_names[prog_name][comm_rank] = []
        prog_names[prog_name][comm_rank].append(value)

    # now that the data is queried and sorted, write it out to the file

    cyclestr = str(simCycle)
    cyclestr = cyclestr.rjust(5,'0')
    filename="performance.metrics." + cyclestr + ".txt"
    flops_out = open(config["outputdir"] + "/" + filename,'w')
    flops_out.write("Process Index, Memory HWM, Total FLOPS, Latest FLOPS, Program Name, Program Index, MPI Rank\n")
    s1 = config["events"]["timer_scaling"][0]
    s2 = config["events"]["timer_scaling"][1]
    index = -1
    prog_index = -1
    current_time = {}
    current_fp_ops = {}
    for prog_name in prog_names:
        prog_index = prog_index + 1
        for comm_rank in prog_names[prog_name]:
            index = index + 1
            current_fp_ops[index] = prog_names[prog_name][comm_rank][1] * s1
            current_time[index] = prog_names[prog_name][comm_rank][2] * s2
            flops_to_date = current_fp_ops[index] / current_time[index]
            if len(last_fp_ops) > 0:
                tmp = (current_fp_ops[index] - last_fp_ops[index]) 
                tmp2 = (current_time[index] - last_time[index])
                if tmp2 > 0.0:
                    # compute flops from lastest timestep
                    flops_in_last_timestep = tmp / tmp2
                else:
                    if last_time[index] > 0.0:
                        # compute flops from previous timestep
                        flops_in_last_timestep = last_fp_ops[index] / last_time[index]
                    else:
                        # something weird is happening...
                        flops_in_last_timestep = 0.0
            else:
                # compute flops from first timestep
                flops_in_last_timestep = flops_to_date
            # write the sorted data to the file
            flops_out.write(str(index) + ", ")
            flops_out.write(str(prog_names[prog_name][comm_rank][0]) + ", ")
            flops_out.write(str(flops_to_date) + ", ")
            flops_out.write(str(flops_in_last_timestep) + ", ")
            flops_out.write(prog_name + ", ")
            flops_out.write(str(prog_index) + ", ")
            flops_out.write(str(comm_rank) + "\n")

    flops_out.close()
    stream_file="performance.metrics.txt"
    stream_out = open(config["outputdir"] + "/" + stream_file,'w')
    stream_out.write(filename + "\n")
    stream_out.close()
    last_time = current_time
    last_fp_ops = current_fp_ops

    return filename
#
#end:def generateVTKFile(...)

def printf(format, *args):
    sys.stdout.write(format % args)

###############################################################################
###############################################################################

if __name__ == "__main__":
    sosScatterplotGenerator()
    #end:FILE


