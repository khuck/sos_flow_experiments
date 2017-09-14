#!/bin/bash -e
set -x

tau=1
impact_makefile=Makefile
reader_script=compile
export COMPILER=ftn 
export ADIOSDIR=${ADIOS_DIR}

if [ ${tau} == 1 ] ; then
    export TAU_ROOT=/ccs/proj/csc143/CODAR_Demo/titan.gnu/tau
    export TAU_ARCH=craycnl
    export TAU_CONFIG=tau-gnu-mpi-pthread-pdt-sos-adios
	export TAU_MAKEFILE=${TAU_ROOT}/${TAU_ARCH}/lib/Makefile.${TAU_CONFIG}
	echo ${TAU_MAKEFILE}
	export PATH=${TAU_ROOT}/${TAU_ARCH}/bin:${PATH}
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
