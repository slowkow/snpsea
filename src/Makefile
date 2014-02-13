# Makefile for SNPsea
# Kamil Slowikowski
# February 12, 2014

# Compiler.
CXX = g++

# Compiler flags.
#DEBUG = -ggdb
CXXFLAGS = $(DEBUG) -w -O3 -m64 -static -fopenmp -std=c++0x
CXXFLAGS += -I ../ext/intervaltree -I ../ext/eigen -DEIGEN_NO_DEBUG
LIB = -static -lz -lpthread -lm -ldl -lgsl

# Header files, source files, objects, binary.
HDR = ezOptionParser.h zfstream.h snpsea.h common.h
SRC = option.cpp zfstream.cpp data.cpp
OBJ = $(SRC:.cpp=.o)
BIN = ../bin/snpsea

# The compiler must be at least 4.6 because we use C++0x features.
GCC_VERSION := $(shell $(CXX) -dumpversion | awk '{print $$1>=4.6?"1":"0"}')

all : $(BIN)

$(BIN) : $(OBJ)

$(OBJ) : $(HDR)

$(BIN) : $(INTERVAL_TREE_PATH) $(EIGEN_PATH)
	mkdir -p ../bin
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJ) $(LIB)

# Make an object file for each C++ source file.
%.o : %.cpp
ifneq "$(GCC_VERSION)" "1"
	$(error You must use GCC 4.6.3 or greater)
else
	$(CXX) $(CXXFLAGS) -c $<
endif

.PHONY: clean

# Create the folder if it doesn't exist.
../ext :
	mkdir $@

# Clone the repositories SNPsea depends on.
../ext/intervaltree :
	git clone 'https://github.com/slowkow/intervaltree.git' $@

../ext/eigen :
	hg clone 'https://bitbucket.org/eigen/eigen' $@

clean:
	rm -f $(OBJ) $(BIN)
