# Library locations
SOS?=$(HOME)/src/sos_flow/install
ADIOS2?=$(HOME)/src/adios2_sst/install_nompi
JSON?=$(HOME)/src/json-3.5.0
TCMALLOC?=$(HOME)/install/gperftools/2.7
JEMALLOC?=$(HOME)/install/jemalloc/5.1.0

# special ops
#SPECIAL=-fsanitize=address -fsanitize=undefined -Wall -Werror
#SPECIAL=-DUSE_TAU_STUBS
#SPECIAL=-L$(JEMALLOC)/lib -ljemalloc 

# Declaration of variables
CC = g++
LINKER = g++
CC_FLAGS = -Wall -Werror -std=c++11 -Ofast -g \
-I$(JSON)/include -I$(SOS)/include \
$(shell $(ADIOS2)/bin/adios2-config --cxx-flags) \
$(SPECIAL)
LD_FLAGS = -Wall -Werror -std=c++11 -Ofast -g \
-L$(SOS)/lib -lsos -Wl,-rpath,$(SOS)/lib \
$(shell $(ADIOS2)/bin/adios2-config --cxx-libs) \
$(SPECIAL)
 
# File names
EXEC = extract
SOURCES = $(wildcard *.cpp) taustubs/tautimer.cpp
HEADERS = $(wildcard *.hpp)
OBJECTS = $(SOURCES:.cpp=.o)
 
# Main target
$(EXEC): $(OBJECTS)
	$(LINKER) $(OBJECTS) -o $(EXEC) $(LD_FLAGS)
	# codesign --force --deep --sign - ./$(EXEC)
 
# To obtain object files
%.o: %.cpp $(HEADERS) Makefile
	$(CC) -c $(CC_FLAGS) $< -o $@
 
# To remove generated files
clean:
	rm -rf $(EXEC) $(OBJECTS) tautrace.bp*

test: $(EXEC)
	time ./$(EXEC)
	# $(ADIOS2)/../build/bin/bpls -A -d ./tautrace.bp
	# $(ADIOS2)/../build/bin/bpls -D -d ./tautrace.bp
