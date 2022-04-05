PREFIX:= /usr/local/

CXX ?= g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl

# htslib/libhts.a nlopt-2.5.0/install/lib/libnlopt.a
# -pg
FLAGS = -Ihtslib -Inloptlib/include -O3 $(LDFLAGS) 

CXXFLAGS := $(FLAGS) $(CXXFLAGS) -std=c++14 

PROGRAMS = DamMet
CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

.PHONY = all clean clean_simple test 

all: $(PROGRAMS)

%.o: %.cpp htslib/libhts.a nloptlib/lib/libnlopt.a
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d

-include $(OBJ:.o=.d)

DamMet: $(OBJ) htslib/libhts.a nloptlib/lib/libnlopt.a
	$(CXX) $^ $(FLAGS) -o DamMet $(LIBS) 


htslib/libhts.a: 
	make -C htslib

nloptlib/lib/libnlopt.a:
	cd nlopt/ && mkdir -p build/  && cd build/ && cmake -DCMAKE_INSTALL_PREFIX=../../nloptlib -DNLOPT_OCTAVE=Off -DNLOPT_MATLAB=Off -DNLOPT_GUILE=Off -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_LIBDIR=lib .. && make && make install && cd ../..

clean:
	rm -f DamMet
	rm -f ./*.o ./*.d 

test: DamMet
	./DamMet
