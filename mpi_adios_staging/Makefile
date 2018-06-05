RM=rm -f
ifeq ($(tau),1)
    CPPFLAGS=-g -Og -I. -std=c++11 $(shell adios_config -c)
    LDLIBS=-lpthread $(shell tau_cxx.sh -tau:showlibs) $(shell adios_config -l)
else
    CPPFLAGS=-g -Og -I. -std=c++11 $(shell adios_config -c)
    LDLIBS=-lpthread $(shell adios_config -l)
endif
LDFLAGS=-g -Og -std=c++11
CC=mpicc
CXX=mpicxx
F90=mpif90

SRCS=$(wildcard *.cpp)
HDRS=$(wildcard *.h)
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))
TARGET=main
NRANKS=4
ITERATIONS=10
WRITE_ITERATIONS=2

all: $(TARGET)

$(OBJS): $(DEPS)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LDLIBS)
	cp $(TARGET) ../.

%.o %.d : %.cpp $(HDRS)

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

test: $(TARGET)
	$(RM) sosd.* profile.* *.bp *.bp.dir
	mpirun -np $(NRANKS) ./$(TARGET) $(ITERATIONS) $(WRITE_ITERATIONS)

debug: $(TARGET)
	$(RM) sosd.* profile.* *.bp *.bp.dir
	mpirun -np $(NRANKS) xterm -e lldb $(shell pwd)/$(TARGET)

tautest: $(TARGET)
	$(RM) sosd.* profile.* *.bp *.bp.dir
	mpirun -np $(NRANKS) tau_exec -T mpi,pthread,adios -ebs ./$(TARGET) $(ITERATIONS) $(WRITE_ITERATIONS)

clean:
	$(RM) $(OBJS) $(TARGET) sosd.* profile.* *.bp *.bp.dir

distclean: clean
	$(RM) *~ $(DEPS)

remake: distclean
	$(MAKE) -j $(TARGET)

ifneq ("$(wildcard $(DEPS))","")
include $(DEPS)
endif
