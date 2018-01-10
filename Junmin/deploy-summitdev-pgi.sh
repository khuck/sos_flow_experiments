#!/bin/bash -e
set -x

origdir=`pwd`

# Where are we going to run the application?  Has to be writable by
# the compute nodes!
rootdir=/lustre/atlas/proj-shared/csc143/`whoami`/testdir-summitdev

# Make sure the path to testdir exists, even though we'll be deleting it.
if [ ! -d ${rootdir} ] ; then
	mkdir -p ${rootdir}
fi

# Where is SOS installed?
sosdir=/ccs/proj/csc143/tau_mona/summitdev.pgi/sos_flow

clean()
{
	# to start fresh, do this.
	rm -rf ${rootdir}
}

setup()
{
	# make the working directory
	if [ ! -d ${rootdir} ] ; then
		mkdir -p ${rootdir}
	fi

	# copy the executables
	cp ImpactTv1betaAdios/xmain ${rootdir}
	cp readerFull/reader2 ${rootdir}

	# copy the run script
	cp summitdev/ImpactT.in ${rootdir}
	cp summitdev/small-pgi.r ${rootdir}

	# copy the SOS daemons
	cp ${sosdir}/bin/sosd ${rootdir}
	cp ${sosdir}/bin/sosd_stop ${rootdir}

}

submit()
{
	cd ${rootdir}
	qsub ./small-pgi.r
	cd ${origdir}
}

clean
setup
#submit
echo "cd ${rootdir} to run the code"

