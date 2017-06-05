#!/bin/bash -e
set -x

origdir=`pwd`
rootdir=${SCRATCH}/mona/sos_flow_experiments/Junmin/testdir

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
	cp cori-MediumRun/ImpactT.in ${rootdir}
	cp cori-MediumRun/cori-small.r ${rootdir}
	#cp tiny/ImpactT.in ${rootdir}
	#cp tiny/cori-small.r ${rootdir}

	# copy the SOS daemons
	cp ${SCRATCH}/mona/sos_flow/build-cori-icc/bin/sosd ${rootdir}
	cp ${SCRATCH}/mona/sos_flow/build-cori-icc/bin/sosd_stop ${rootdir}
}

submit()
{
	cd ${rootdir}
	sbatch ./cori-small.r
	cd ${origdir}
}

clean
setup
#submit
