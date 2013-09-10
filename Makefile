# Makefile for SNPspec
# 
# Supported platforms: Unix / Linux

# Directory of the target
OUTPUT = snpspec

# Compiler
CXX = g++

# EIGEN library
EIGEN_PATH = /home/kamil/src/eigen-eigen-ffa86ffb5570

# Compiler flags
CXXFLAGS = -w -O3 -m64 -static -fopenmp -std=c++0x -I $(EIGEN_PATH) -DEIGEN_NO_DEBUG
LIB += -static -lz -lpthread -lm -ldl -lgsl

HDR += ezOptionParser.h IntervalTree.h zfstream.h snpspec.h common.h
SRC = option.cpp zfstream.cpp data.cpp
	   
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :
	$(CXX) $(CXXFLAGS) -o $(OUTPUT) $(OBJ) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)

FORCE:

clean: 
	rm -f *.o
