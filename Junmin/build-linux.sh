#!/bin/bash -e
set -x

tau=1
impact_makefile=Makefile
reader_script=compile

if [ ${tau} == 1 ] ; then
	export TAU_CONFIG=-papi-mpi-pthread-sos-adios
	export TAU_ARCH=x86_64
	export TAU_ROOT=$HOME/src/tau2
	export ADIOS_ROOT=$HOME/src/ADIOS/ADIOS-gcc
	export ADIOSDIR=$ADIOS_ROOT
	export TAU_MAKEFILE=${TAU_ROOT}/${TAU_ARCH}/lib/Makefile.tau${TAU_CONFIG}
	export PATH=${TAU_ROOT}/${TAU_ARCH}/bin:${ADIOS_ROOT}/bin:${PATH}
	impact_makefile=Makefile.tau
	reader_script=compile.tau
fi

impact()
{
	cd ImpactTv1betaAdios
	make clean
	make -f ${impact_makefile}
	cd ..
}

reader()
{
	cd readerFull
	rm -f *.o reader2
	./${reader_script}
	cd ..
}

impact
reader
