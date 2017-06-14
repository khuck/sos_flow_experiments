include Makefile.inc
RM=rm -f
CPPFLAGS=-g -I. -I$(SOSDIR)/include -pthread $(CXX_STD_FLAG) $(POWER_FLAG)
LDFLAGS=-g -lpthread $(CXX_STD_FLAG)
LDLIBS=-L$(SOSDIR)/lib -lsos

SRCS=$(wildcard *.cpp)
HDRS=$(wildcard *.h)
OBJS=$(subst .cpp,.o,$(SRCS))
DEPS=$(subst .cpp,.d,$(SRCS))
TARGET=main

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

clean:
	$(RM) $(OBJS) $(TARGET) sosd.*

distclean: clean
	$(RM) *~ $(DEPS)

remake: distclean
	$(MAKE) -j $(TARGET)

ifneq ("$(wildcard $(DEPS))","")
include $(DEPS)
endif
