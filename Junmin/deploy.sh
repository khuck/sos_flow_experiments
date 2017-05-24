#!/bin/bash -e
set -x

origdir=`pwd`
rootdir=/home/khuck/src/MONA/testdir

clean()
{
	# to start fresh, do this.
	rm -rf ${rootdir}
}

setup()
{
	# make the working directory
	if [ ! -d ${rootdir} ] ; then
		mkdir ${rootdir}
	fi

	# copy the executables
	cp ImpactTv1betaAdios/xmain ${rootdir}
	cp readerFull/reader2 ${rootdir}

	# copy the run script
	#cp tiny/ImpactT.in ${rootdir}
	#cp tiny/small.r ${rootdir}

	# copy the SOS daemons
	cp $HOME/src/sos_flow/build-linux/bin/sosd ${rootdir}
	cp $HOME/src/sos_flow/build-linux/bin/sosd_stop ${rootdir}

	# copy the ADIOS library
	cp $HOME/src/chaos/linux-gcc/lib/libenet.so* ${rootdir}
}

submit()
{
	cd ${rootdir}
	#qsub ./small.r
	cd ${origdir}
}

#clean
setup
#submit
