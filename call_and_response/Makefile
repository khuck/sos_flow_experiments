include Makefile.inc
CC=mpicc
CXX=mpicxx
RM=rm -f
CPPFLAGS=-g -I$(SOSDIR)/include -pthread
LDFLAGS=-g -lpthread
LDLIBS=-L$(SOSDIR)/lib -lsos

SRCS=$(wildcard *.cpp)
OBJS=$(subst .cpp,.o,$(SRCS))
TARGET=main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LDLIBS) 

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CC) -M $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

clean:
	$(RM) $(OBJS) $(TARGET) sosd.*

distclean: clean
	$(RM) *~ .depend

