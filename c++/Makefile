.PHONY: clean all info

TARGETS := main
SOURCES := main.cpp mesh.cpp integrals.cpp
OBJS    := $(SOURCES:.cpp=.o) $(SOURCES:.cpp=_d.o) $(SOURCES:.cpp=_p.o)


DEBUGFLAGS:= -g3
RELEASEFLAGS:= -O0 -DNDEBUG
PARALLELFLAGS:= -D_GLIBCXX_PARALLEL -fopenmp -pthread
CXXFLAGS := -Wall -std=c++17
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



integrals.o: integrals.cpp $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(DEBUGFLAGS) integrals.cpp
mesh.o: mesh.cpp integrals.o $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(DEBUGFLAGS) mesh.cpp
main.o: main.cpp mesh.o integrals.o $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(DEBUGFLAGS) main.cpp
TrustRegion.o: TrustRegion.cpp mesh.o $(EXAMPLE_DEPS)
	 @$(CXX) -c $(CXXFLAGS) $(DEBUGFLAGS) TrustRegion.cpp

main: main.o mesh.o integrals.o TrustRegion.o
	 @$(CXX) -o main $(LDFLAGS) $(CXXFLAGS) main.o mesh.o integrals.o TrustRegion.o $(LIBS)
