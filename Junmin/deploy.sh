#!/bin/bash -e
set -x

origdir=`pwd`
# rootdir=/home/khuck/src/MONA/testdir
# rootdir=${origdir}/testdir
rootdir=/lustre/atlas/proj-shared/csc143/khuck/testdir
sosdir=/ccs/proj/csc143/tau/sos_flow

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
	# cp tiny/ImpactT.in ${rootdir}
	# cp tiny/small.r ${rootdir}

	# copy the SOS daemons
	cp ${sosdir}/bin/sosd ${rootdir}
	cp ${sosdir}/bin/sosd_stop ${rootdir}

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
