HOST=$(shell hostname)
ifeq ($(HOST),ktau)
  SOSDIR=$(HOME)/src/sos_flow/build-linux
  CC=mpicc
  CXX=mpicxx
  POWER_FLAG=
else ifeq ($(HOST),ln1)
  SOSDIR=$(HOME)/src/sos_flow/build-linux
  CC=mpicc
  CXX=mpicxx
  POWER_FLAG=
else ifeq ($(NERSC_HOST),cori)
  SOSDIR=$(SCRATCH)/mona/sos_flow/build-cori-${compiler}
  POWER_FLAG=-DAPEX_HAVE_CRAY_POWER
  CC=cc
  CXX=CC
else ifeq ($(BUILDDIR),build-falcon)
  SOSDIR=$(ADIOS_ROOT)
  POWER_FLAG=
  CC=mpicc
  CXX=mpicxx
else
  SOSDIR=$(HOME)/src/sos_flow/build-titan
  CC=cc
  CXX=CC
  POWER_FLAG=
endif
