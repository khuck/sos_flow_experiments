#!/bin/bash -e
set -x

tau=1
impact_makefile=Makefile
reader_script=compile

if [ ${tau} == 1 ] ; then
    if [ ${PE_ENV} == "INTEL" ] ; then
	  export TAU_CONFIG=-intel-mpi-pthread-pdt
    else
	  export TAU_CONFIG=-gnu-mpi-pthread-pdt
    fi
	export TAU_ARCH=craycnl
	export TAU_ROOT=${SCRATCH}/mona/tau2
	export TAU_MAKEFILE=${TAU_ROOT}/${TAU_ARCH}/lib/Makefile.tau${TAU_CONFIG}
	export PATH=${TAU_ROOT}/${TAU_ARCH}/bin:${PATH}
	impact_makefile=Makefile.tau
	reader_script=compile.tau
fi

impact()
{
	cd ImpactTv1betaAdios
	#make clean
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

if [ ${PE_ENV} == "INTEL" ] ; then
  export COMPILER="ftn"
else
  export COMPILER="ftn -ffree-line-length-0"
fi
export ADIOSDIR=${SCRATCH}/mona/adios/git-${compiler}
impact
export COMPILER="cc"
reader
