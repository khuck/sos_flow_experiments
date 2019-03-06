# Library locations
SOS=$(HOME)/src/sos_flow/install
ADIOS2=$(HOME)/src/adios2_greg/install_nompi
JSON=$(HOME)/src/json-3.5.0

# Declaration of variables
CC = g++
CC_FLAGS = -w -I$(JSON)/include -I$(SOS)/include \
$(shell $(ADIOS2)/bin/adios2-config --cxx-flags) -std=c++11 -O0 -g
LD_FLAGS = -L$(SOS)/lib -lsos -Wl,-rpath,$(SOS)/lib \
$(shell $(ADIOS2)/bin/adios2-config --cxx-libs)
 
# File names
EXEC = extract
SOURCES = $(wildcard *.cpp)
HEADERS = $(wildcard *.hpp)
OBJECTS = $(SOURCES:.cpp=.o)
 
# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) -o $(EXEC) $(LD_FLAGS)
	codesign --force --deep --sign - ./$(EXEC)

 
# To obtain object files
%.o: %.cpp $(HEADERS) Makefile
	$(CC) -c $(CC_FLAGS) $< -o $@
 
# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)

test: $(EXEC)
	./$(EXEC)
