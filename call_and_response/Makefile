include Makefile.inc
CC=mpicc
CXX=mpicxx
RM=rm -f
CPPFLAGS=-g -I$(SOSDIR)/include
LDFLAGS=-g
LDLIBS=-L$(SOSDIR)/lib -lsos

SRCS=$(wildcard *.cpp)
OBJS=$(subst .cpp,.o,$(SRCS))
TARGET=main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LDLIBS) 

depend: .depend

.depend: $(SRCS)
	$(RM) ./.depend
	$(CXX) $(CPPFLAGS) -MM $^>>./.depend;

clean:
	$(RM) $(OBJS) $(TARGET) sosd.*

distclean: clean
	$(RM) *~ .depend

include .depend
