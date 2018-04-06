RM=rm -f
CPPFLAGS=-g -Og -I. -std=c++11
LDFLAGS=-g -Og -std=c++11
LDLIBS=-lpthread
CC=mpicc
CXX=mpicxx
F90=mpif90

SRCS=$(wildcard *.cpp)
HDRS=$(wildcard *.h)
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))
TARGET=main
NRANKS=36

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
	mpirun -np $(NRANKS) ./$(TARGET)

debug: $(TARGET)
	mpirun -np $(NRANKS) xterm -e lldb $(shell pwd)/$(TARGET)

tautest: $(TARGET)
	mpirun -np $(NRANKS) tau_exec -T mpi -ebs ./$(TARGET)

clean:
	$(RM) $(OBJS) $(TARGET) sosd.* profile.*

distclean: clean
	$(RM) *~ $(DEPS)

remake: distclean
	$(MAKE) -j $(TARGET)

ifneq ("$(wildcard $(DEPS))","")
include $(DEPS)
endif
