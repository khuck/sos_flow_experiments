# Library locations
SOS=$(HOME)/src/sos_flow/install
ADIOS2=$(HOME)/src/adios2_sst/install_nompi
JSON=$(HOME)/src/json-3.5.0

# Declaration of variables
CC = g++
#CC_FLAGS = -w -I$(SOS)/include $(shell $ADIOS2/bin/adios2-config --cxx-flags)
CC_FLAGS = -w -I$(JSON)/include -I$(SOS)/include -I$(ADIOS2)/include -std=c++11
 
# File names
EXEC = extract
SOURCES = $(wildcard *.cpp)
OBJECTS = $(SOURCES:.cpp=.o)
 
# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC)
 
# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@
 
# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)