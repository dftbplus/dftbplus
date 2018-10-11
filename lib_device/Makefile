# Source, Includes, Library, Output Defines
CXXFLAGS	= -I/path.to.magma/include -I/path.to.cuda/include
SRC_F  		= magma_device_info.f90
SRC_C  		= functions.cpp
OBJ_F  		= $(SRC_F:.f90=.o)
OBJ_C  		= $(SRC_C:.cpp=.o)
LDFLAGS   	= -L/path.to.magma/lib -lmagma
TARGET   	= libCountDevice.a

# Compiler Defines
CXX    = g++
FC     = gfortran
FFLAGS = -O3 -fopenmp

# Compile F90 Source File into Object File
%.o: %.f90
		$(FC) -c $*.f90

# Compile CXX Source File into Object File
%.o: %.cpp
		$(CXX) -c $(CXXFLAGS) $(LDFLAGS) $*.cpp

# Generate the staic library
$(TARGET): $(OBJ_F) $(OBJ_C)
		ar rcs $(TARGET) $(OBJ_F) $(OBJ_C)

# Clean Up Objects, Module, and Library
clean:
		$(RM) $(TARGET) $(OBJ_F) $(OBJ_C) *.mod

