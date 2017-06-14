HOST=$(shell hostname)
CXX_STD_FLAG=-std=c++11 
ifeq ($(HOST),ktau)
  SOSDIR=$(HOME)/src/sos_flow/build-linux
  CC=mpicc
  CXX=mpicxx
  POWER_FLAG=
else ifeq ($(HOST),ln1)
  SOSDIR=$(HOME)/src/sos_flow/build-linux
  CC=mpicc
  CXX=mpicxx
  POWER_FLAG=-DAPEX_HAVE_POWERCAP_POWER
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
  CXX_STD_FLAG=-std=c++0x 
else
  SOSDIR=$(HOME)/src/sos_flow/build-titan
  CC=cc
  CXX=CC
  POWER_FLAG=
endif
