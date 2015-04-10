CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wno-deprecated-declarations -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 # -fopenmp
INC = -I/usr/local/include -I$(HOME)/local/include
CXXFLAGS += $(INC)
COMPILE.c = $(CXX) $(CXXFLAGS)
LDFLAGS = -lz -lbz2 -pthread -lpthread #-fopenmp
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


utils.o: utils.cpp utils.hpp 
bam_parser.o: bam_parser.cpp reference_info_builder.hpp ribomap_profiler.hpp bam_parser.hpp utils.hpp
reference_info_builder.o: reference_info_builder.cpp reference_info_builder.hpp utils.hpp
ribomap_profiler.o: ribomap_profiler.cpp reference_info_builder.hpp ribomap_profiler.hpp bam_parser.hpp
footprint_generator.o: footprint_generator.cpp reference_info_builder.hpp ribomap_profiler.hpp transcript_model.hpp utils.hpp lattice.hpp
lattice.o: lattice.cpp lattice.hpp utils.hpp math_utils.hpp
