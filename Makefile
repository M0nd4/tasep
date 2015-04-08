CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations
COMPILE.c = $(CXX) $(CXXFLAGS)
LDFLAGS = 

all: build/tasep_playground
clean:
	rm -f build/*.o build/tasep_playground

build/tasep_playground: build/tasep_playground.o build/lattice.o
	$(COMPILE.c) $(LDFLAGS) -o build/tasep_playground build/tasep_playground.o build/lattice.o
	chmod u+x $@

build/tasep_playground.o: src/tasep_playground.cpp src/lattice.hpp src/utils.hpp
	$(COMPILE.c) -c -o build/tasep_playground.o src/tasep_playground.cpp

build/lattice.o: src/lattice.cpp src/lattice.hpp src/utils.hpp src/math_utils.hpp
	$(COMPILE.c) -c -o build/lattice.o src/lattice.cpp
