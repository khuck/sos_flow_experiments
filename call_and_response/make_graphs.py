#!/usr/bin/env python
import os
import sys
import sqlite3
import numpy as np
import pylab as pl
import time
import signal
import random

conn = None
min_timestamp = None
graphs = [None, None, None, None, None, None]
axises = [None, None, None, None, None, None]

def open_connection(sqlite_file):
    global conn
    # check for file to exist
    print ("Checking for file: ", sqlite_file)
    while not os.path.exists(sqlite_file):
        print ("Waiting on file: ", sqlite_file)
        time.sleep(1)

    print("Connecting to: ", sqlite_file)
    # Connecting to the database file
    #conn = sqlite3.connect(sqlite_file)
    #fd = os.open(sqlite_file, os.O_RDONLY)
    #conn = sqlite3.connect('/dev/fd/%d' % fd)
    #url = 'file:' + sqlite_file + '?mode=ro'
    #url = 'file:' + sqlite_file
    #conn = sqlite3.connect(url, uri=True)
    conn = sqlite3.connect(sqlite_file)
    conn.isolation_level=None
    c = conn.cursor()
    #c.execute('PRAGMA journal_mode=WAL;')
    #c.execute('PRAGMA synchronous   = ON;')
    #c.execute('PRAGMA cache_size    = 31250;')
    #c.execute('PRAGMA cache_spill   = FALSE;')
    #c.execute('PRAGMA temp_store    = MEMORY;')
    return c

def signal_handler(signal, frame):
    print("Detected ctrl-C...exiting.")
    print("Closing connection to database.")
    conn.close()
    sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

def try_execute(c, statement, parameters=None):
    try:
        if parameters:
            c.execute(statement,parameters);
        else:
            c.execute(statement);
    except sqlite3.Error as e:
        print("database error...", e.args[0])

def make_index(c):
    sql_statement = ("create index foo2 on tblvals(guid);")
    #print("Executing query")
    try_execute(c,sql_statement);

def get_ranks(c):
    sql_statement = ("select distinct comm_rank, process_id from tblpubs where title not like 'system monitor' and title not like 'process monitor: %' order by comm_rank;")
    #print("Executing query")
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    ranks = np.array([x[0] for x in all_rows])
    procs = np.array([x[1] for x in all_rows])
    ranklen = len(ranks)
    if ranklen > 10:
        smallranks = [0]
        for i in range(1,4):
            candidate = random.randrange(1, ranklen-1)
            while candidate in smallranks:
                candidate = random.randrange(1, ranklen-1)
            smallranks.append(candidate)
        smallranks.append(int(ranklen-1))
        smallranks2 = []
        smallprocs2 = []
        for index in smallranks:
            smallranks2.append(ranks[index])
            smallprocs2.append(procs[index])
        return np.array(sorted(smallranks2)), np.array(sorted(smallprocs2))
    else:
        return ranks, procs

def get_nodes(c):
    sql_statement = ("select distinct node_id, min(comm_rank) from tblpubs group by node_id order by comm_rank;")
    #print("Executing query")
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    nodes = np.array([x[0] for x in all_rows])
    ranks = np.array([x[1] for x in all_rows])
    nodelen = len(nodes)
    if nodelen > 10:
        smallnodes = [0]
        for i in range(1,4):
            candidate = random.randrange(1, nodelen-1)
            while candidate in smallnodes:
                candidate = random.randrange(1, nodelen-1)
            smallnodes.append(candidate)
        smallnodes.append(int(nodelen-1))
        smallnodes2 = []
        smallranks2 = []
        for index in sorted(smallnodes):
            smallnodes2.append(nodes[index])
            smallranks2.append(ranks[index])
        return np.array(smallnodes2), np.array(smallranks2)
    else:
        return nodes,ranks

def get_min_timestamp(c):
    global min_timestamp
    sql_statement = ("select min(time_pack) from tblvals;")
    #print("Executing query")
    try_execute(c,sql_statement);
    all_rows = c.fetchall()
    ts = np.array([x[0] for x in all_rows])
    min_timestamp = ts[0]
    print("min timestamp: ", min_timestamp)

def do_chart(subplot, c, ranks, ranks2, group_column, metric, plot_title, y_label, graph, axes):
    global min_timestamp
    newplot = False
    if not graph:
        newplot = True
        graph = {}
    index = 0
    for r in ranks:
        sql_statement = ("SELECT distinct tbldata.name, tblvals.val, tblvals.time_pack, tblpubs.comm_rank FROM tblvals INNER JOIN tbldata ON tblvals.guid = tbldata.guid INNER JOIN tblpubs ON tblpubs.guid = tbldata.pub_guid WHERE tblvals.guid IN (SELECT guid FROM tbldata WHERE tbldata.name LIKE '" + metric + "') AND tblpubs." + group_column)
        """
        if isinstance(r, int):
            sql_statement = (sql_statement + " = " + str(r) + " order by tblvals.time_pack;")
        else:
            sql_statement = (sql_statement + " like '" + r + "' order by tblvals.time_pack;")
        """
        if ranks2 == []:
            sql_statement = (sql_statement + " = " + str(r) + " and tblvals.val > 0 order by tblvals.time_pack;")
        else:
            sql_statement = (sql_statement + " like '" + r + "' and tblvals.val > 0 order by tblvals.time_pack;")

        #params = [metric,r]
        print "Executing query: ", sql_statement,
        try_execute(c, sql_statement)
        print "Done. "

        #print("Fetching rows.")
        all_rows = c.fetchall()
        if len(all_rows) <= 0:
            print("Error: query returned no rows.",)
            print(sql_statement, params)

        #print("Making numpy array of: metric_values")
        metric_values = np.array([max(x[1],0) for x in all_rows])
        #print("Making numpy array of: pack_time")
        pack_time = np.array([x[2]-min_timestamp for x in all_rows])

        #print("len(pack_time) == ", len(pack_time))
        #print("len(metric_values) == ", len(metric_values))

        #print("Plotting: x=pack_time, y=metric_values")
        if newplot:
            axes = pl.subplot(subplot)
            axes.set_title(plot_title);
            graph[r] = (pl.plot(pack_time, metric_values, marker='*', linestyle='-', label=str(r))[0])
            axes.set_autoscale_on(True) # enable autoscale
            axes.autoscale_view(True,True,True)
            pl.legend(prop={'size':6})
            pl.ylabel(y_label)
            pl.xlabel("Timestamp")
        else:
            graph[r].set_data(pack_time, metric_values)
            axes.relim()        # Recalculate limits
            axes.autoscale_view(True,True,True) #Autoscale
        index = index + 1
    return graph,axes

def do_derived_chart(subplot, c, ranks, group_column, metric1, metric2, plot_title, y_label, graph, axes):
    global min_timestamp
    newplot = False
    if not graph:
        newplot = True
        graph = {}
    for r in ranks:
        sql_statement = ("SELECT tbldata.name, cast(tblvals.val as float), tblvals.time_pack FROM tblvals INNER JOIN tbldata ON tblvals.guid = tbldata.guid INNER JOIN tblpubs ON tblpubs.guid = tbldata.pub_guid WHERE tblvals.guid IN (SELECT guid FROM tbldata WHERE tbldata.name LIKE '" + metric1 + "') AND tblpubs." + group_column + " = " + str(r) + " order by tblvals.time_pack;")

        # print("Executing query 1")
        params = [metric1,r]
        try_execute(c,sql_statement);
        # print("Fetching rows.")
        all_rows1 = c.fetchall()
        if len(all_rows1) <= 0:
            print("Error: query returned no rows.",)
            print(sql_statement, params)

        sql_statement = ("SELECT tbldata.name, cast(tblvals.val as float), tblvals.time_pack FROM tblvals INNER JOIN tbldata ON tblvals.guid = tbldata.guid INNER JOIN tblpubs ON tblpubs.guid = tbldata.pub_guid WHERE tblvals.guid IN (SELECT guid FROM tbldata WHERE tbldata.name LIKE '" + metric2 + "') AND tblpubs." + group_column + " = " + str(r) + " order by tblvals.time_pack LIMIT " + str(len(all_rows1)) + ";")
        # print("Executing query 2")
        params = [metric2,r]
        try_execute(c,sql_statement);
        # print("Fetching rows.")
        all_rows2 = c.fetchall()
        if len(all_rows2) <= 0:
            print("Error: query returned no rows.",)
            print(sql_statement, params)

        # print("Making numpy array of: metric_values")
        metric_values = np.array([x[1]/y[1] for x,y in zip(all_rows1,all_rows2)])
        # print(metric_values)
        # print("Making numpy array of: pack_time")
        pack_time = np.array([x[2]-min_timestamp for x in all_rows1])
        # print(pack_time)
        if len(metric_values) > len(pack_time):
            np.resize(metric_values, len(pack_time))
        elif len(pack_time) > len(metric_values):
            np.resize(pack_time, len(metric_values))

        # print("len(pack_time) == ", len(pack_time))
        # print("len(metric_values) == ", len(metric_values))

        # print("Plotting: x=pack_time, y=metric_values")
        if newplot:
            axes = pl.subplot(subplot)
            axes.set_title(plot_title);
            graph[r] = (pl.plot(pack_time, metric_values, marker='*', linestyle='-', label=str(r))[0])
            axes.set_autoscale_on(True) # enable autoscale
            axes.autoscale_view(True,True,True)
            pl.legend(prop={'size':6})
            pl.ylabel(y_label)
            pl.xlabel("Timestamp")
        else:
            graph[r].set_data(pack_time, metric_values)
            axes.relim()        # Recalculate limits
            axes.autoscale_view(True,True,True) #Autoscale
        pl.draw()
    return graph,axes

def docharts(c,nodes,noderanks,ranks,procs):
    # rows, columns, figure number for subplot value
    graphs[0],axises[0] = do_chart(321, c, ranks, [], "comm_rank", "Iteration","Time per iteration","Seconds", graphs[0], axises[0])
    graphs[1],axises[1] = do_chart(322, c, nodes, noderanks, "node_id", "CPU System%","CPU System","CPU Utilization (%)", graphs[1], axises[1])
    #graph3,axes3 = do_chart(323, c, ranks, "comm_rank", "status:VmRSS%","Mean memory footprint (KB)","Kilobytes", graphs[0], axises[0])
    #graphs[2],axises[2] = do_chart(323, c, procs, [], "process_id", "status:VmRSS%","Resident Set Size","Kilobytes", graphs[2], axises[2])
    graphs[2],axises[2] = do_chart(323, c, procs, [], "process_id", "Matrix Size","Data Set Size","Cells", graphs[2], axises[2])
    graphs[3],axises[3] = do_chart(324, c, nodes, noderanks, "node_id", "CPU User%","CPU User","CPU Utilization (%)", graphs[3], axises[3])
    graphs[4],axises[4] = do_chart(325, c, procs, [], "process_id", "status:VmHWM","High Water Mark","Kilobytes", graphs[4], axises[4])
    graphs[5],axises[5] = do_chart(326, c, nodes, noderanks, "node_id", "Package-0 Energy","Package-0 Energy","Package-0 Energy", graphs[5], axises[5])
    #graph6,axes6 = do_derived_chart(326, c, ranks, "comm_rank", "%TAU::0::exclusive_TIME::MPI_Waitall()%","%TAU::0::calls::MPI_Waitall()%","MPI_Waitall() ","MPI_Waitall()", [], [])

def main(arguments):
    global min_timestamp
    # name of the sqlite database file
    sqlite_file = arguments[0]

    # open the connection
    c = open_connection(sqlite_file)

    # get the number of ranks
    make_index(c)
    ranks,procs = get_ranks(c)
    while ranks.size == 0:
        time.sleep(1)
        ranks,procs = get_ranks(c)
    print ("ranks: ", ranks)
    # get the number of nodes
    nodes,noderanks = get_nodes(c)
    while nodes.size == 0:
        time.sleep(1)
        nodes,noderanks = get_nodes(c)
    print ("nodes: ", nodes)
    get_min_timestamp(c)
    #resize the figure
    # Get current size
    fig_size = pl.rcParams["figure.figsize"]
    # Set figure width to 12 and height to 9
    fig_size[0] = 12
    fig_size[1] = 9
    pl.rcParams["figure.figsize"] = fig_size
    pl.ion()
    docharts(c,nodes,noderanks,ranks,procs)
    print("Closing connection to database.")
    # Closing the connection to the database file
    conn.close()
    pl.tight_layout()
    while True:
        pl.pause(30.0)
        print("Updating chart...")
        # open the connection
        c = open_connection(sqlite_file)
        docharts(c,nodes,noderanks,ranks,procs)
        print("Closing connection to database.")
        # Closing the connection to the database file
        conn.close()

    print("Done.")

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
