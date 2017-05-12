CC=mpicc
CXX=mpicxx
RM=rm -f
CPPFLAGS=-g -I/home/khuck/src/sos_flow/build-linux/include
LDFLAGS=-g
LDLIBS=-L/home/khuck/src/sos_flow/build-linux/lib -lsos

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
