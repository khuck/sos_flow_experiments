#!/usr/bin/env python
import os
import sys
import sqlite3
import numpy as np
import yaml
from collections import  OrderedDict

conn = None
comm_dict = {}

##############################################################################
# A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise
# Martin Ester, Hans-Peter Kriegel, Jorg Sander, Xiaowei Xu
# dbscan: density based spatial clustering of applications with noise

import math

UNCLASSIFIED = False
NOISE = None

def _dist(p,q):
    return math.sqrt(np.power(p-q,2).sum())

def _eps_neighborhood(p,q,eps):
    return _dist(p,q) < eps

def _region_query(m, point_id, eps):
    n_points = m.shape[1]
    seeds = []
    for i in range(0, n_points):
        if _eps_neighborhood(m[:,point_id], m[:,i], eps):
            seeds.append(i)
    return seeds

def _expand_cluster(m, classifications, point_id, cluster_id, eps, min_points):
    seeds = _region_query(m, point_id, eps)
    if len(seeds) < min_points:
        classifications[point_id] = NOISE
        return False
    else:
        classifications[point_id] = cluster_id
        for seed_id in seeds:
            classifications[seed_id] = cluster_id
            
        while len(seeds) > 0:
            current_point = seeds[0]
            results = _region_query(m, current_point, eps)
            if len(results) >= min_points:
                for i in range(0, len(results)):
                    result_point = results[i]
                    if classifications[result_point] == UNCLASSIFIED or \
                       classifications[result_point] == NOISE:
                        if classifications[result_point] == UNCLASSIFIED:
                            seeds.append(result_point)
                        classifications[result_point] = cluster_id
            seeds = seeds[1:]
        return True
        
def dbscan(m, eps, min_points):
    """Implementation of Density Based Spatial Clustering of Applications with Noise
    See https://en.wikipedia.org/wiki/DBSCAN
    
    scikit-learn probably has a better implementation
    
    Uses Euclidean Distance as the measure
    
    Inputs:
    m - A matrix whose columns are feature vectors
    eps - Maximum distance two points can be to be regionally related
    min_points - The minimum number of points to make a cluster
    
    Outputs:
    An array with either a cluster id number or dbscan.NOISE (None) for each
    column vector in m.
    """
    cluster_id = 1
    n_points = m.shape[1]
    classifications = [UNCLASSIFIED] * n_points
    for point_id in range(0, n_points):
        point = m[:,point_id]
        if classifications[point_id] == UNCLASSIFIED:
            if _expand_cluster(m, classifications, point_id, cluster_id, eps, min_points):
                cluster_id = cluster_id + 1
    return classifications

def test_dbscan():
    m = np.matrix('1 1.2 0.8 3.7 3.9 3.6 10; 1.1 0.8 1 4 3.9 4.1 10')
    eps = 0.5
    min_points = 2
    assert dbscan(m, eps, min_points) == [1, 1, 1, 2, 2, 2, None]

def dbscan_vector(vector, eps, min_points):
    m = np.matrix(vector)
    """
    # get the distances between all points
    distances = []
    for i in range(0,len(vector)):
        for j in range(i+1,len(vector)):
            distances.append(np.abs(vector[i] - vector[j]))
    # make an autobinning histogram of the values
    """
    print vector
    hist,binedges = np.histogram(vector,bins='rice')
    print hist
    print binedges
    right = 0.0
    left = 0.0
    foundright = False
    for i in xrange(len(hist)-1,-1,-1):
        if hist[i] == 0:
            right = binedges[i+1]
            foundright = True
            min_points = hist[i+1]
        elif hist[i] > 0 and foundright:
            left = binedges[i+1]
            eps = right - left
            break
    print("eps: ",eps)
    clusters = dbscan(m, eps, min_points)
    # found out which cluster has largest values
    subclust = {}
    maxval = 0.0
    maxc = 0
    for c,v in zip(clusters,vector):
        if c != None:
            if c not in subclust:
                subclust[c] = []
            subclust[c].append(v)
            if v > maxval:
                maxval = v
                maxc = c
                #print maxval,maxc
    return clusters,len(subclust),maxc,np.min(subclust[maxc])
##############################################################################

# Make a connection to the SQLite3 database
def open_connection(filename):
    global conn
    # check for file to exist
    print ("Checking for file: ", sqlite_file)
    while not os.path.exists(sqlite_file):
        print ("Database file does not exist: ", sqlite_file)
        sys.exit(1)

    print("Connecting to: ", sqlite_file)
    # Connecting to the database file
    conn = sqlite3.connect(sqlite_file)
    #fd = os.open(sqlite_file, os.O_RDONLY)
    #conn = sqlite3.connect('/dev/fd/%d' % fd)
    url = 'file:' + sqlite_file + '?mode=ro'
    #url = 'file:' + sqlite_file
    #conn = sqlite3.connect(url, uri=True)
    conn.isolation_level=None
    conn.text_factory = bytes
    c = conn.cursor()
    return c

# wrapper around queries, for error handling
def try_execute(c, statement, parameters=None):
    success = False
    print(statement)
    while not success:
        try:
            if parameters:
                c.execute(statement,parameters);
            else:
                c.execute(statement);
            success = True
            break;
        except sqlite3.Error as e:
            #print("database error...", e.args[0])
            success = False

# find all of the ranks participating in the simulation, sending data over SOS
def get_ranks(c,application):
    sql_statement = ("select distinct guid,comm_rank from tblPubs where prog_name like '%" + application + "%' order by comm_rank;")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    pub_guids = np.array([x[0] for x in all_rows])
    ranks = np.array([x[1] for x in all_rows])
    ranklen = len(ranks)
    #print ("pub_guids: ", pub_guids)
    #print ("ranks: ", ranks)
    return pub_guids, ranks

# Get the GUIDs for the MPI Collective events, for this publisher GUID (i.e. for this rank).
def get_mpi_collective_guid(pg):
    sql_statement = ("select guid from tblData where pub_guid = " + str(pg) + " and name like 'MPI collective%';")
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    guid = np.array([x[0] for x in all_rows])
    # print "MPI collective exchange guid: ", guid
    tmpstr = ''
    for g in guid:
        tmpstr = tmpstr + str(g) + ','
    tmpstr = tmpstr[:-1]
    return tmpstr

# Get the GUIDs for the ADIOS WRITE  events, for this publisher GUID (i.e. for this rank).
def get_adios_write_guid(pg):
    sql_statement = ("select guid from tblData where pub_guid = " + str(pg) + " and name like 'ADIOS %';")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    guid = np.array([x[0] for x in all_rows])
    # print "MPI collective exchange guid: ", guid
    tmpstr = ''
    for g in guid:
        tmpstr = tmpstr + str(g) + ','
    tmpstr = tmpstr[:-1]
    return tmpstr

# get the time range of interest.  It is three of the last four iterations, ignoring
# the last one that could be bogus.
def get_advance_step_range(application, limit, order):
    sql_statement = "select time_pack from viewCombined where value_name like 'TAU_EVENT::adios_advance_step%' and comm_rank = 0 and prog_name like '" + application + "' order by time_pack " + order + " limit " + str(limit) + ";"
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    times = np.array([x[0] for x in all_rows])
    return(times[0],times[limit-1])

# get the time range of interest.  It is three of the last four iterations, ignoring
# the last one that could be bogus.
def get_time_range(pg):
    limit = 2
    sql_statement = ("select distinct tbldata.name, tblvals.val, tblvals.time_pack from tblvals left outer join tbldata on tblvals.guid = tbldata.guid where tbldata.name like 'ADIOS close%' and tbldata.pub_guid = " + str(pg) + " order by time_pack desc limit " + str(limit) + ";")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    closes = np.array([x[2] for x in all_rows])
    #print("closes:",closes)
    sql_statement = ("select distinct tbldata.name, tblvals.val, tblvals.time_pack from tblvals left outer join tbldata on tblvals.guid = tbldata.guid where tbldata.name like 'ADIOS open%' and tbldata.pub_guid = " + str(pg) + " order by time_pack desc limit " + str(limit) + ";")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    opens = np.array([x[2] for x in all_rows])
    #print("opens:",opens)
    #print("range:",str(closes[limit-1]),str(opens[0]))
    return(closes[limit-1],opens[0])

# get the time range of interest for ADIOS calls.  It is the first iteration.
def get_time_range_adios(pg):
    limit = 1
    sql_statement = ("select distinct tbldata.name, tblvals.val, tblvals.time_pack from tblvals left outer join tbldata on tblvals.guid = tbldata.guid where tbldata.name like 'ADIOS close%' and tbldata.pub_guid = " + str(pg) + " order by time_pack asc limit " + str(limit) + ";")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    closes = np.array([x[2] for x in all_rows])
    #print("closes:",closes)
    sql_statement = ("select distinct tbldata.name, tblvals.val, tblvals.time_pack from tblvals left outer join tbldata on tblvals.guid = tbldata.guid where tbldata.name like 'ADIOS open%' and tbldata.pub_guid = " + str(pg) + " order by time_pack asc limit " + str(limit) + ";")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    opens = np.array([x[2] for x in all_rows])
    #print("opens:",opens)
    #print("range:",str(closes[limit-1]),str(opens[0]))
    return(opens[0],closes[0])

# Get the last n start and stop timestamps for all MPI collective events for this rank. 
# First, sort them by "most recent first", then take those n rows and reverse the order 
# so we have them in chronological order.
def get_mpi_exchanges(mpi_guid,limit,time_start,time_end):
    sql_statement = ("select start, end, name from (select time_pack-val as start, time_pack as end, name from tblVals inner join tbldata on tblvals.guid = tbldata.guid where tbldata.guid in (" + str(mpi_guid) + ") and time_pack > " + repr(time_start) + " and time_pack < " + repr(time_end) + " order by end DESC limit " + limit + ") order by end ASC;")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    starts = np.array([x[0] for x in all_rows])
    ends = np.array([x[1] for x in all_rows])
    names = np.array([x[2] for x in all_rows])
    return starts,ends,names

# Get the last n start and stop timestamps for all MPI collective events for this rank. 
# First, sort them by "most recent first", then take those n rows and reverse the order 
# so we have them in chronological order.
def get_trace_events(application, time_start,time_end):
    start_cond = ""
    end_cond = ""
    if time_start != 0:
        start_cond = " and time_pack > " + repr(time_start)
    if time_end != 0:
        end_cond = " and time_pack <= " + repr(time_end)
    sql_statement = ("select value_name, time_pack-value, time_pack from viewCombined where value_name like 'TAU_EVENT::%' and comm_rank = 0 and prog_name like '%" + application + "'" + start_cond + end_cond + " order by time_pack asc;")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    names = np.array([x[0] for x in all_rows])
    starts = np.array([x[1] for x in all_rows])
    ends = np.array([x[2] for x in all_rows])
    return starts,ends,names

# Get the last n start and stop timestamps for all ADIOS write events for this rank. 
# First, sort them by "most recent first", then take those n rows and reverse the order 
# so we have them in chronological order.
def get_adios_writes(mpi_guid,limit,time_start,time_end):
    sql_statement = ("select start, end, name from (select time_pack-val as start, time_pack as end, name from tblVals inner join tbldata on tblvals.guid = tbldata.guid where tbldata.guid in (" + str(mpi_guid) + ") and time_pack > " + repr(time_start) + " and time_pack < " + repr(time_end) + " order by end DESC limit " + limit + ") order by end ASC;")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    starts = np.array([x[0] for x in all_rows])
    ends = np.array([x[1] for x in all_rows])
    names = np.array([x[2] for x in all_rows])
    return names

# Just take the difference between the two arrays, start[n] - end[n-1] for all n.
# The resulting array should have n-1 elements.
def find_gap_between_timestamps(starts,ends):
    if len(ends) == 0 or len(starts) == 0:
        return [0];
    durations = np.zeros(len(ends)-1)
    for i in range(1,len(ends)):
        durations[i-1] = starts[i]-ends[i-1]
    return durations

# For validation.
def find_nth_percentile(durations,n):
    print n, "th percentile:", np.percentile(durations,n)

# For validation.
def reject_outliers(data, m=2):
    means = data[abs(data - np.mean(data)) < m * np.std(data)]
    medians = data[abs(data - np.median(data)) < m * np.std(data)]
    #print np.min(means), np.mean(means), np.median(means), np.max(means)
    #print np.min(medians), np.mean(medians), np.median(medians), np.max(medians)
    return means,medians

def get_first_arrival(data):
    return np.min(data)

def get_last_arrival(data):
    return np.max(data)

def setup_yaml():
    """ https://stackoverflow.com/a/8661021 """
    represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)    

def get_comm_yaml(application,op_list):
    global comm_dict
    sql_statement = ("select distinct tbldata.name from tblvals left outer join tbldata on tblvals.guid = tbldata.guid left outer join tblpubs on tbldata.pub_guid = tblpubs.guid where (tbldata.name like 'MPI_Comm%' or tbldata.name like 'MPI_Cart%') and tblpubs.prog_name like '%" + application + "';")
    #print(sql_statement)
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    commands = np.array([x[0] for x in all_rows])
    for com in commands:
        #print(com)
        tmp = com.split(' ',1)
        name = tmp[0]
        tmp2 = tmp[1].rsplit(' ',1)
        args = tmp2[0]
        result = tmp2[1]
        my_dict = OrderedDict()
        my_dict["type"] = str(name)
        tmp = args.strip('(').split(',',1)
        comm_value = tmp[0].strip(')')
        comm_name = "MPI_COMM_WORLD";
        # do we have MPI_COMM_WORLD yet?
        if len(comm_dict) > 0:
            comm_name = "MPI_COMM_" + str(len(comm_dict))
        if comm_value not in comm_dict:
            comm_dict[str(comm_value)] = comm_name
        else:
            comm_name = comm_dict[comm_value]
        #print(comm_value, comm_name)
        my_dict["args"] = str(args.replace(comm_value,comm_name))
        my_dict["name"] = str(comm_name)
        # handle the result
        comm_value = result
        comm_name = "MPI_COMM_" + str(len(comm_dict))
        if comm_value not in comm_dict:
            comm_dict[str(comm_value)] = comm_name
        else:
            comm_name = comm_dict[comm_value]
        #print(comm_value, comm_name)
        #my_dict["result_value"] = str(comm_value)
        my_dict["result"] = str(comm_name)
        op_list.append(my_dict)
    if len(commands) == 0:
        sql_statement = ("select distinct value_name from viewCombined where value_name like 'TAU_EVENT::MPI%' and prog_name like '%" + application + "' order by time_pack asc limit 1;")
        try_execute(c,sql_statement);
        all_rows = c.fetchall()
        commands = np.array([x[0] for x in all_rows])
        for com in commands:
            tokens = com.split(" ")
            # get the last token, that's the communicator
            comm_value = tokens[-1]
            comm_dict[comm_value] = "MPI_COMM_WORLD"

def get_adios_yaml(names,op_list):
    short = ''
    my_dict = OrderedDict()
    my_dict["type"] = "adios phase"
    my_dict["bp_file"] = "pookie_extraction.bp"
    my_dict["comm"] = "MPI_COMM_WORLD"
    my_vars = []
    for n in names:
        print n
        var_dict = OrderedDict()
        tmp = n.split('(',1)
        args = tmp[1].split(')',1)
        print args
        arg = args[0].split(',',3)
        if len(arg) > 2:
            var_dict["name"] = str(arg[0])
            var_dict["type"] = str(arg[1])
            ndims = int(arg[2])
            rest = arg[3].rsplit(',',1)
            #print(arg[3],rest)
            if ndims > 0:
                if rest[0].startswith('[') and rest[0].endswith(']'):
                    rest[0] = rest[0][1:-1]
                #var_dict["local_dims"] = str(rest[0].strip('[',1).strip(']',1).strip(','))
                var_dict["local_dims"] = str(rest[0].strip(','))
            else:
                var_dict["local_dims"] = ""
                var_dict["value"] = int(rest[1])
            my_vars.append(var_dict)
    my_dict["vars"] = my_vars
    op_list.append(my_dict)

def get_combined_yaml(starts,ends,names,op_list,maxv):
    global comm_dict
    # anything within 50% of the biggest cluster is a real compute region
    threshold = float(maxv) * 0.5
    previous_end = 0
    for s,e,n in zip(starts,ends,names):
        if "MPI" in n:
            #print n
            if previous_end > 0:
                possible_compute = float(s) - float(previous_end)
                if possible_compute > threshold:
                    # do some compute
                    my_dict = OrderedDict()
                    my_dict["type"] = str("compute")
                    my_dict["time"] = str(possible_compute)
                    op_list.append(my_dict)
            previous_end = e
            short = ''
            my_dict = OrderedDict()
            if 'collective exchange' in n:
                short = n.replace('TAU_EVENT::MPI collective exchange ','MPI_')
                tmp = short.split(' ',2)
                my_dict["type"] = str(tmp[0])
                my_dict["num_bytes"] = int(tmp[1].strip('(').rstrip(')'))
                comm_name = str(tmp[2].strip())
                if str(tmp[2].strip()) in comm_dict:
                    comm_name = comm_dict[str(tmp[2].strip())]
                my_dict["comm"] = str(comm_name)
            if 'collective synchronize' in n:
                short = n.replace('TAU_EVENT::MPI collective synchronize ','MPI_')
                tmp = short.split(' ',1)
                my_dict["type"] = str(tmp[0])
                comm_name = "MPI_COMM_WORLD"
                if str(tmp[1].strip()) in comm_dict:
                    comm_name = comm_dict[str(tmp[1].strip())]
                my_dict["comm"] = str(comm_name)
            op_list.append(my_dict)
        else:
            var_dict = OrderedDict()
            tokens = n.split('(',1)
            var_dict["type"] = tokens[0].replace("TAU_EVENT::","")
            args = tokens[1].rstrip(')')
            if len(args) > 0:
                if "comm: " in args:
                    for key in comm_dict:
                        args = args.replace(str(key), comm_dict[key])
                tokens = args.split(', ')
                for t in tokens:
                    t = t.strip()
                    left_right = t.split(':',1)
                    var_dict[left_right[0].strip()] = left_right[1].strip().replace("'","")
                op_list.append(var_dict)
     
def get_mpi_yaml(starts,ends,names,op_list,maxv):
    global comm_dict
    # anything within 50% of the biggest cluster is a real compute region
    threshold = float(maxv) * 0.5
    previous_end = 0
    for s,e,n in zip(starts,ends,names):
        #print n
        if previous_end > 0:
            possible_compute = float(s) - float(previous_end)
            if possible_compute > threshold:
                # do some compute
                my_dict = OrderedDict()
                my_dict["type"] = str("compute")
                my_dict["time"] = str(possible_compute)
                op_list.append(my_dict)
        previous_end = e
        short = ''
        my_dict = OrderedDict()
        if 'collective exchange' in n:
            short = n.replace('MPI collective exchange ','MPI_')
            tmp = short.split(' ',2)
            my_dict["type"] = str(tmp[0])
            my_dict["num_bytes"] = int(tmp[1].strip('(').rstrip(')'))
            comm_name = str(tmp[2].strip())
            comm_name = "MPI_COMM_WORLD"
            if str(tmp[2].strip()) in comm_dict:
                comm_name = comm_dict[str(tmp[2].strip())]
            my_dict["comm"] = str(comm_name)
        if 'collective synchronize' in n:
            short = n.replace('MPI collective synchronize ','MPI_')
            tmp = short.split(' ',1)
            my_dict["type"] = str(tmp[0])
            comm_name = "MPI_COMM_WORLD"
            if str(tmp[1].strip()) in comm_dict:
                comm_name = comm_dict[str(tmp[1].strip())]
            my_dict["comm"] = str(comm_name)
        op_list.append(my_dict)

setup_yaml()

if len(sys.argv) != 3:
    print("Usage: brute-force.py <database> <executable-name>")
    print("Example: ./brute-force.py sosd.00000.db reader2")
    sys.exit(0)

# name of the sqlite database file
sqlite_file = sys.argv[1]
application = sys.argv[2]

# open the connection
c = open_connection(sqlite_file)

# First, get the start time of the first and second advance_step calls
first, second = get_advance_step_range(application,2, "asc")

# Then get the list of events that happened before the first advance step call
init_starts,init_ends,init_names = get_trace_events(application,0,first)

# Then get the list of events that happened between the first and second advance_step calls
loop_starts,loop_ends,loop_names = get_trace_events(application,first,second)

# remove the overlap to get the initialization phase
init_len = len(init_names)
loop_len = len(loop_names)
new_init_len = init_len - loop_len
init_starts = init_starts[:new_init_len]
init_ends = init_ends[:new_init_len]
init_names = init_names[:new_init_len]

# get the last advance_step call
last, next_to_last = get_advance_step_range(application,2, "desc")

# get the finalization steps
fini_starts,fini_ends,fini_names = get_trace_events(application,last,0)

# get the publisher guids, ranks
pub_guids,ranks = get_ranks(c, application)

# populate the YAML
master_dict=OrderedDict()
master_dict["ranks"] = len(ranks)
init_op_list = []
loop_op_list = []
fini_op_list = []
get_comm_yaml(application,init_op_list)

# How long are the compute windows?
durations = find_gap_between_timestamps(loop_starts,loop_ends)
#print(maxv=durations[0]
if len(durations) > 1:
    clusters,numclust,maxc,maxv = dbscan_vector(durations, 0.5, 2)
print numclust
#print durations
print clusters,numclust,maxc,maxv

# add the init phase to YAML
get_combined_yaml(init_starts,init_ends,init_names,init_op_list,maxv)
master_dict["init_op_list"] = init_op_list
# add the loop phase to YAML

get_combined_yaml(loop_starts,loop_ends,loop_names,loop_op_list,maxv)
master_dict["loop_op_list"] = loop_op_list
# add the fini phase to YAML
get_combined_yaml(fini_starts,fini_ends,fini_names,fini_op_list,maxv)
master_dict["fini_op_list"] = fini_op_list

# write the YAML
outfile = open('output.yaml','w')
yaml.dump(master_dict,outfile,default_flow_style=True)
outfile.close

