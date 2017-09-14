# sos_flow_experiments

This repository contains one example that uses SOS and TAU for runtime
monitoring of a parallel applications and workflows (a pipeline, really).

# Configuration

TAU, SOS, etc. are installed on titan in /ccs/proj/csc143/tau.  To build
and run this example on titan, get this code from github and do the following:

1) run build-titan.sh to build xmain and reader2, the two applications in the
pipeline.

2) modify deploy-titan.sh to change the location of ${rootdir}.  Make sure that
it is a writable directory from the compute nodes!

3) run deploy-titan.sh to copy everything to the execution directory (${rootdir})

4) cd to the execution directory

5) qsub small.r

# Output

There should be the following relevant output:

* both.out - the output from the small.r script
* xmain.out - the output from xmain
* read2.out - the output from reader2
* profiles_xmain/ - a directory containing the TAU profiles from xmain
* profiles_reader2/ - a directory containing the TAU profiles from reader2
* sosd.00000.db - a sqlite3 database containing the aggregated performance data
* sosd.0000[1-5].db - a sqlite3 database containing the performance data from 
each listener

# Todo: 

We will soon be adding an SOS analysis node that will extract the performance
data at runtime, and generate VTK files for visualization.
