#!/usr/bin/env python
#

import sys
import time
import pprint as pp
from ssos import SSOS

def demonstrateSOS():
    SOS = SSOS()

    print "Initializing SOS..."
    SOS.init()

    sql_string = "SELECT * FROM tblVals WHERE rowid < 4000;"
    print "Sending this query to the SOS daemon: "
    print "    " + sql_string
    results, col_names = SOS.query(sql_string)
    print "Results:"
    print "    Output rows....: " + str(len(results))    #pp.pprint(results)
    print "    Column count...: " + str(len(col_names)) 
    print "    Column names...: "# + str(col_names)    #pp.pprint(col_names)
    pp.pprint(col_names)
    print ""
    print "Finalizing..."
    SOS.finalize();
    print "   ...DONE!"
    print 

if __name__ == "__main__":
    demonstrateSOS()




