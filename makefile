CXX = g++
CXXFLAGS = -std=c++11
OBJECTS = BoxTreeNode.o FMMUtilities.o NbodyFMM.o
EXE = NbodyFMM

all : $(EXE)

$(EXE) : $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJECTS):
	$(CXX) $(CXXFLAGS) -c $<

NbodyFMM.o: NbodyFMM.cpp BoxTreeNode.hpp FMMUtilities.hpp
BoxTreeNode.o: BoxTreeNode.cpp BoxTreeNode.hpp
FMMUtilities.o: FMMUtilities.cpp FMMUtilities.hpp

clean:
	rm -f *.o *.out

.PHONY: all clean