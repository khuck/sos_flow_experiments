#!/bin/bash -e
set -x

tau=1
impact_makefile=Makefile
reader_script=compile

if [ ${tau} == 1 ] ; then
	export TAU_CONFIG=-mpi-pthread-pdt
	export TAU_ARCH=x86_64
	export TAU_ROOT=$HOME/src/tau2
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

impact
reader
