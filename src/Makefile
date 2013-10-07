# Makefile for SNPspec
# 
# Supported platforms: Unix / Linux

# Compiler
CXX = g++

EIGEN_PATH = ../eigen-eigen-ffa86ffb5570
INTERVAL_TREE_PATH = ../intervaltree

INCLUDE = -I $(EIGEN_PATH) -I $(INTERVAL_TREE_PATH)

DEBUG = #-ggdb

# Compiler flags
CXXFLAGS = $(DEBUG) -w -O3 -m64 -static -fopenmp -std=c++0x 
CXXFLAGS += $(INCLUDE) -DEIGEN_NO_DEBUG
LIB = -static -lz -lpthread -lm -ldl -lgsl

HDR = ezOptionParser.h zfstream.h snpspec.h common.h
SRC = option.cpp zfstream.cpp data.cpp
OBJ = $(SRC:.cpp=.o)

OUTPUT = ../bin/snpspec

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR) $(INTERVALTREE)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean: 
	rm -f $(OBJ) $(OUTPUT)
