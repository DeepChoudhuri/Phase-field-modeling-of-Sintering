CXX=/usr/bin/g++
#CXX=/usr/bin/gcc
#CXX=clang

CXXFLAGS=-c -Wall -O2 -std=c++11
#CXXFLAGS=-c -O2 -std=c++0x -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=1 -std=c++11
FFTWLIB=-lfftw3

EXE = ceq.x

all: $(EXE)

ceq.x: main.o PFUtilities.o CoupledEqFDSinter.o Parameters.o InitialStructures.o
	$(CXX) $(FFTWLIB) -O2 main.o PFUtilities.o CoupledEqFDSinter.o Parameters.o InitialStructures.o -o ceq.x

main.o: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp

PFUtilities.o: PFUtilities.cpp PFUtilities.h
	$(CXX) $(CXXFLAGS) PFUtilities.cpp

CoupledEqSinter.o: CoupledEqFDSinter.cpp CoupledEqFDSinter.h
	$(CXX) $(CXXFLAGS) CoupledEqFDSinter.cpp

Parameters.o: Parameters.cpp Parameters.h
	$(CXX) $(CXXFLAGS) Parameters.cpp

InitialStructures.o: InitialStructures.cpp InitialStructures.h
		$(CXX) $(CXXFLAGS) InitialStructures.cpp


run:
	time ./$(EXE) > log

clean:
	rm -f *o *.x *.vtk log* *.txt
	clear

vtkclean:
	rm -f *.vtk
	clear

clear:
	rm -f *o *.x *.vtk log* *.txt
	clear
