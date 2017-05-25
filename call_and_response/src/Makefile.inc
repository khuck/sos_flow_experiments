HOST=$(shell hostname)
ifeq ($(HOST),ktau)
  SOSDIR=$(HOME)/src/sos_flow/build-linux
  CC=mpicc
  CXX=mpicxx
else ifeq ($(HOST),ln1)
  SOSDIR=$(HOME)/src/sos_flow/build-linux
  CC=mpicc
  CXX=mpicxx
else ifeq ($(NERSC_HOST),cori)
  SOSDIR=$(SCRATCH)/mona/sos_flow/build-cori-${compiler}
  CC=cc
  CXX=CC
else
  SOSDIR=$(HOME)/src/sos_flow/build-titan
  CC=cc
  CXX=CC
endif
