include Makefile.inc
RM=rm -f
CPPFLAGS=-g -I$(SOSDIR)/include -pthread -std=c++11
LDFLAGS=-g -lpthread -std=c++11
LDLIBS=-L$(SOSDIR)/lib -lsos

SRCS=$(wildcard *.cpp)
HDRS=$(wildcard *.h)
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))
TARGET=main

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LDLIBS) 

%.o %.d : %.cpp $(HDRS)

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

clean:
	$(RM) $(OBJS) $(TARGET) sosd.*

distclean: clean
	$(RM) *~ $(DEPS)

ifneq ("$(wildcard $(DEPS))","")
include $(DEPS)
endif
