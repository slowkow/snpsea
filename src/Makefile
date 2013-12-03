# Makefile for SNPsea
# Supported platforms: Unix / Linux

# Compiler.
CXX = g++

# Paths to external resources.
EIGEN_PATH = ../ext/eigen
EIGEN_URL = https://bitbucket.org/eigen/eigen

INTERVAL_TREE_PATH = ../ext/intervaltree
INTERVAL_TREE_URL = https://github.com/slowkow/intervaltree.git

# Compiler flags.
#DEBUG = -ggdb
INCLUDE = -I $(EIGEN_PATH) -I $(INTERVAL_TREE_PATH)
CXXFLAGS = $(DEBUG) -w -O3 -m64 -static -fopenmp -std=c++0x
CXXFLAGS += $(INCLUDE) -DEIGEN_NO_DEBUG
LIB = -static -lz -lpthread -lm -ldl -lgsl

# Header files, source files, objects, binary.
HDR = ezOptionParser.h zfstream.h snpspec.h common.h
SRC = option.cpp zfstream.cpp data.cpp
OBJ = $(SRC:.cpp=.o)
BIN = ../bin/snpspec

all : $(BIN)

$(BIN) : intervaltree eigen
	$(CXX) $(CXXFLAGS) -o $(BIN) $(OBJ) $(LIB)

$(OBJ) : $(HDR)

# Make an object file for each C++ source file.
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

$(BIN) : $(OBJ)

.PHONY: intervaltree eigen clean

intervaltree: $(INTERVAL_TREE_PATH)
	test -d $< || git clone $(INTERVAL_TREE_URL) $<

eigen: $(EIGEN_PATH)
	test -d $< || hg clone $(EIGEN_URL) $<

clean:
	rm -f $(OBJ) $(BIN)
