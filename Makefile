OPENMPFLAG ?= -fopenmp
CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 ${OPENMPFLAG}
INC = -I/usr/local/include -I$(HOME)/local/include -I.
CXXFLAGS += $(INC)
COMPILE.c = $(CXX) $(CXXFLAGS)
NVCC = nvcc
LDFLAGS = -lz -lbz2 -pthread -lpthread ${OPENMPFLAG}
OUTPUT_OPTION = -o $@

all: footprint_generator tasep_playground
clean:
	rm -f *.o footprint_generator tasep_playground

footprint_generator: footprint_generator.o reference_info_builder.o ribomap_profiler.o transcript_model.o utils.o lattice.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@

tasep_playground: lattice.o tasep_playground.o utils.o
	$(COMPILE.c) $(OUTPUT_OPTION) $^ $(LDFLAGS)
	chmod u+x $@



ribo: lattice_ribo.o utils.o tasep_playground.o
	$(COMPILE.c) -o tasep_playground lattice_ribo.o utils.o tasep_playground.o $(LDFLAGS) -lcudart -L $(CUDA_HOME)/lib64 

lattice_ribo.o: lattice_ribo.cu lattice.hpp
	$(NVCC) -c -arch=sm_20 -o lattice_ribo.o lattice_ribo.cu



single: lattice_single.o utils.o tasep_playground.o
	$(COMPILE.c) -o tasep_playground lattice_single.o utils.o tasep_playground.o $(LDFLAGS)

lattice_single.o: lattice_single.cpp lattice.hpp
	$(COMPILE.c) -c -o lattice_single.o lattice_single.cpp



approx_cpu: lattice_approx.o utils.o tasep_playground.o
	$(COMPILE.c) -o tasep_playground lattice_approx.o utils.o tasep_playground.o $(LDFLAGS)

lattice_approx.o: lattice_approx.cpp lattice.hpp
	$(COMPILE.c) -c -o lattice_approx.o lattice_approx.cpp



tasep_playground.o: tasep_playground.cpp
	$(COMPILE.c) -c -o tasep_playground.o tasep_playground.cpp

utils.o: utils.cpp utils.hpp 
	$(COMPILE.c) -c -o utils.o utils.cpp


bam_parser.o: bam_parser.cpp reference_info_builder.hpp ribomap_profiler.hpp bam_parser.hpp utils.hpp
reference_info_builder.o: reference_info_builder.cpp reference_info_builder.hpp utils.hpp
ribomap_profiler.o: ribomap_profiler.cpp reference_info_builder.hpp ribomap_profiler.hpp bam_parser.hpp
footprint_generator.o: footprint_generator.cpp reference_info_builder.hpp ribomap_profiler.hpp transcript_model.hpp utils.hpp lattice.hpp
lattice.o: lattice.cpp lattice.hpp utils.hpp math_utils.hpp
