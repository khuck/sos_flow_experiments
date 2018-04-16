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

# Global attr, representing the "previous" frame
previous_attr = dict()
cached_results_dim = None
cached_selected_fields = None
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

# Iterate over the running aggregators, and query them.
def queryAllAggregators(sql):
    global sosPort
    global aggregators
    global SOS
    allResults = None
    # Iterate over the aggregators, running the query against each one.
    # This could be done in parallel...
    for a in aggregators:
        print a,sosPort,sql,
        sys.stdout.flush()
        results, col_names = SOS.query(sql, a, sosPort)
        print "\n...done."
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
        os.environ[str(var)] = config["sosd"][var]
    if "periodic" not in config:
        config["periodic"] = False

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
    global cached_selected_fields
    global cachedName
    global config
    global last_time
    global last_fp_ops
    #####
    #
    # Get the list of field names we will use to build a custom query.
    #
    # NOTE: To filter out SOS_VAL_TYPE_STRING fields, add in:
    #            ... += "WHERE value_type != 3"
    if cached_selected_fields == None:
        sqlFieldNames = "SELECT DISTINCT name FROM tbldata where val_type != 3 order by name;"
        results, col_names = queryAllAggregators(sqlFieldNames)
        tmplist = [el[0] for el in results]
        # we queried against multiple servers, so get the superset of field names
        tmpset = set (tmplist)
        # change it back to a list
        selectedFields = list(tmpset)
        cached_selected_fields = selectedFields
    else:
        selectedFields = cached_selected_fields
    name_count = len(selectedFields)

    printf("(%d fields)\n", name_count)

    #
    # NOTE: Debug output...
    #
    #print "Selected " + str(name_count) + " unique names:"
    #for name in selectedFields:
    #    print "    " + str(name)
    #print ""
    #
    #####


    # NOTE: VisIt will give errors when you are projecting a field that is
    #       not present in all of its data sets. If that were going to be
    #       a problem, i.e. for production code, it might be good to filter
    #       the list of fields to only those present in ALL data sets, and
    #       provide the user with a list of fields that were not included
    #       (for this reason) ... or include fields but with a special
    #       value for 'not present'.
    
    # Get the PMI data
    start = time.time()
    totalstart = start

    # get the max frame for each publisher that finished by that time
    sqlFieldNames = "drop table if exists lastframe;"
    results, col_names = queryAllAggregators(sqlFieldNames)
    sqlFieldNames = "create temporary table lastframe (pg unsigned big int, g unsigned big int, f integer);"
    results, col_names = queryAllAggregators(sqlFieldNames)
    if config["periodic"]:
        sqlFieldNames = "insert into lastframe select pub_guid, tbldata.guid, " + str(simCycle) + " as mframe from tbldata where tbldata.name = '" + cachedName + "' group by pub_guid order by pub_guid;"
    else:
        if config["aggregators"]["runtime"]:
            operator = ">="
        else:
            operator = "="
        sqlFieldNames = "insert into lastframe select pub_guid, tbldata.guid, max(frame) as mframe from tbldata left outer join tblvals on tbldata.guid = tblvals.guid where tbldata.name = '" + cachedName + "' and frame " + operator + " " + str(simCycle) + " group by pub_guid order by pub_guid;"
    results, col_names = queryAllAggregators(sqlFieldNames)
    end = time.time()
    print (end-start), "seconds for max(frame) query"

    # Get the frame-specific data...

    # do the memory first

    start = time.time()
    # Make sure we sort by prog_name, comm_rank!
    sqlValsToColByRank = "select pg, val, f from lastframe left outer join tbldata on lastframe.pg = tbldata.pub_guid and tbldata.name = '" + str(config["events"]["counters"][0]) + "' left outer join tblvals on tbldata.guid = tblvals.guid and tblvals.frame = lastframe.f left outer join tblpubs on pg = tblpubs.guid order by tblpubs.prog_name, tblpubs.comm_rank;"
    results, col_names = queryAllAggregators(sqlValsToColByRank)
    print "query done"
    end = time.time()
    print (end-start), "seconds for frame query"

    cyclestr = str(simCycle)
    cyclestr = cyclestr.rjust(5,'0')

    # filename="performance.memory." + cyclestr + ".txt"
    # memory_out = open(config["outputdir"] + "/" + filename,'w')
    # memory_out.write("\"Rank\", \"" + str(config["events"]["counters"][0]) + "\"\n")
    index = 0
    memory_data = {}
    for r in results:
        # memory_out.write(str(index) + ", " + str(r[1]) + "\n")
        memory_data[index] = r[1]
        index = index + 1
    # memory_out.close()
    # bpfile="performance.memory.streaming.bp"
    # bp_out = open(config["outputdir"] + "/" + bpfile,'w')
    # bp_out.write(filename + "\n")
    # bp_out.close()

    # do the timer PAPI metric next

    start = time.time()
    # Make sure we sort by prog_name, comm_rank!
    sqlValsToColByRank = "select pg, val, f from lastframe left outer join tbldata on lastframe.pg = tbldata.pub_guid and tbldata.name = '" + str(config["events"]["timers"][0]) + "' left outer join tblvals on tbldata.guid = tblvals.guid and tblvals.frame = lastframe.f left outer join tblpubs on pg = tblpubs.guid order by tblpubs.prog_name, tblpubs.comm_rank;"
    dividends, col_names = queryAllAggregators(sqlValsToColByRank)
    print "query done"
    end = time.time()
    print (end-start), "seconds for frame query"

    # do the timer time metric next

    start = time.time()
    # Make sure we sort by prog_name, comm_rank!
    sqlValsToColByRank = "select pg, val, f, prog_name, comm_rank from lastframe left outer join tbldata on lastframe.pg = tbldata.pub_guid and tbldata.name = '" + str(config["events"]["timers"][1]) + "' left outer join tblvals on tbldata.guid = tblvals.guid and tblvals.frame = lastframe.f left outer join tblpubs on pg = tblpubs.guid order by tblpubs.prog_name, tblpubs.comm_rank;"
    divisors, col_names = queryAllAggregators(sqlValsToColByRank)
    print "query done"
    end = time.time()
    print (end-start), "seconds for frame query"

    filename="performance.metrics." + cyclestr + ".txt"
    flops_out = open(config["outputdir"] + "/" + filename,'w')
    flops_out.write("Index, Memory HWM, Total FLOPS, Latest FLOPS, Program Name, Program Index, MPI Rank\n")
    s1 = config["events"]["timer_scaling"][0]
    s2 = config["events"]["timer_scaling"][1]
    index = 0
    prog_index = 0
    prog_name = str(divisors[0][3])
    current_time = {}
    current_fp_ops = {}
    for d1,d2 in zip(dividends,divisors):
        current_fp_ops[index] = float(d1[1]) * s1
        current_time[index] = float(d2[1]) * s2
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
        flops_out.write(str(index) + ", ")
        flops_out.write(str(memory_data[index]) + ", ")
        flops_out.write(str(flops_to_date) + ", ")
        flops_out.write(str(flops_in_last_timestep) + ", ")
        flops_out.write(str(d2[3]) + ", ")
        if str(d2[3]) != prog_name:
            prog_index = prog_index + 1
            prog_name = str(d2[3])
        flops_out.write(str(prog_index) + ", ")
        flops_out.write(str(d2[4]) + "\n")
        index = index + 1
    flops_out.close()
    bpfile="performance.metrics.txt"
    bp_out = open(config["outputdir"] + "/" + bpfile,'w')
    bp_out.write(filename + "\n")
    bp_out.close()
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


