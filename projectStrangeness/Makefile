# Makefile for compiling a C++ program

# Compiler
CXX = g++

# Compiler flags
CXXFLAGS = -std=c++11 -Wall -Wextra `root-config --cflags`

# Source files
SRCS = ssbar_generate.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

# Include options
INC = -I$(PYTHIA8)/include 

#LDFLAGS
LDFLAGS=`root-config --libs`

# Linking
LINK = -L$(PYTHIA8)/lib -lpythia8

# Executable name
EXEC = ssbar_generate

# Default target
all: $(EXEC)

# Rule for compiling object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INC) -c $< -o $@

# Rule for linking object files into executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LINK) $(OBJS) -o $(EXEC)

# Clean rule
clean:
	$(RM) $(OBJS) $(EXEC)

