.PHONY: clean all info

TARGETS := main
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


helpers.o: helpers.cpp $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS) helpers.cpp
integrals.o: integrals.cpp helpers.o $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS) integrals.cpp
mesh.o: mesh.cpp helpers.o integrals.o $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS) mesh.cpp
main.o: main.cpp mesh.o integrals.o helpers.o $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS) main.o main.cpp


main: main.o mesh.o integrals.o helpers.o
	 @$(CXX) -o main $(LDFLAGS) $(CXXFLAGS) main.o mesh.o integrals.o helpers.o  $(LIBS)
