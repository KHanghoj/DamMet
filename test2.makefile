PREFIX:= /usr/local/

CXX ?= g++

LIBS = -lz -lm -lbz2 -llzma -lpthread -lcurl





# htslib/libhts.a nlopt-2.5.0/install/lib/libnlopt.a
FLAGS = -Ihtslib -Inlopt-2.5.0/install/include -O3 $(LDFLAGS)

CXXFLAGS := $(FLAGS) $(CXXFLAGS)

PACKAGE_VERSION = 1.0.2

PROGRAMS = DamMet
CXXSRC = $(wildcard *.cpp)
OBJ = $(CXXSRC:.cpp=.o)

all: $(PROGRAMS)

%.o: %.cpp
	$(CXX) -c  $(CXXFLAGS) $*.cpp
	$(CXX) -MM $(CXXFLAGS) $*.cpp >$*.d

version.h:
	echo '#define DamMet_VERSION "$(PACKAGE_VERSION)"' > $@

-include $(OBJ:.o=.d)

DamMet: version.h $(OBJ) htslib/libhts.a nlopt-2.5.0/install/lib/libnlopt.a
	$(CXX) $(FLAGS) -o DamMet $(LIBS) $^

htslib/libhts.a:
	git clone https://github.com/samtools/htslib.git
	cd htslib  && make && cd ../

nlopt-2.5.0/install/lib/libnlopt.a: v2.5.0.tar.gz
	cd nlopt-2.5.0/ && mkdir -p build/  && cd build/ && cmake -DCMAKE_INSTALL_PREFIX=../install -DNLOPT_OCTAVE=Off -DNLOPT_MATLAB=Off -DNLOPT_GUILE=Off -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_LIBDIR=../install/lib .. && make && make install && cd ../..

v2.5.0.tar.gz:
	wget https://github.com/stevengj/nlopt/archive/v2.5.0.tar.gz
	tar xf v2.5.0.tar.gz
