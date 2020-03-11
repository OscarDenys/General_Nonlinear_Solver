.PHONY: clean all info

TARGETS := main main_d main mlmc mlmc_d
SOURCES := main.cpp mesh.cpp integrals.cpp helpers.cpp
OBJS    := $(SOURCES:.cpp=.o) $(SOURCES:.cpp=_d.o) $(SOURCES:.cpp=_p.o)


DEBUGFLAGS:= -g
RELEASEFLAGS:= -O0 -DNDEBUG
PARALLELFLAGS:= -D_GLIBCXX_PARALLEL -fopenmp -pthread
CXXFLAGS := -Wall -std=c++1z
LDFLAGS  :=
LIBS =   -lstdc++
OPENMPLIB =-L/usr/lib -fopenmp -pthread

EXAMPLE_DEPS = Makefile


all: $(TARGETS)

clean:
	/bin/rm -f $(OBJS) $(TARGETS)

info:
	@echo Compiler:	 CXX	  = $(CXX)
	@echo Compile command: COMPILE.cc  = $(COMPILE.cc)
	@echo Link command:    LINK.cc     = $(LINK.cc)

mesh.o: mesh.cpp $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS) mesh.o mesh.cpp
integrals.o: integrals.cpp mehs.o $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS)  -o integrals.o integrals.cpp mesh.o
main.o: main.cpp mesh.cpp integrals.cpp $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS)  -o main.o main.cpp


main: main.o mesh.o integrals.o
	 @$(CXX) -o main $(LDFLAGS) $(CXXFLAGS) main.o mesh.o integrals.o  $(LIBS)
