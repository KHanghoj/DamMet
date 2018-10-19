

DamMet: deammeth.cc myargparser.cc nlopt-2.5.0/install/include/nlopt.hpp htslib/libhts.a
	g++  -O3  deammeth.cc myargparser.cc ./htslib/libhts.a nlopt-2.5.0/install/lib/libnlopt.a -std=c++11 -I nlopt-2.5.0/install/include  -I./htslib/ -lz -lm -lbz2 -llzma -o DamMet -lpthread
	# g++  -O3  deammeth.cc myargparser.cc ./htslib/libhts.a -std=c++11 -I./htslib/ -lz -lm -lbz2 -llzma -lnlopt_cxx -o DamMet -lpthread

nlopt-2.5.0/install/include/nlopt.hpp:
	wget https://github.com/stevengj/nlopt/archive/v2.5.0.tar.gz
	tar xvf v2.5.0.tar.gz
	cd nlopt-2.5.0/ && mkdir -p build/  && cd build/ && cmake -DCMAKE_INSTALL_PREFIX=../install -DBUILD_SHARED_LIBS=OFF .. && make && make install && cd ../..

htslib/libhts.a:
	git clone https://github.com/samtools/htslib.git
	cd htslib && make -j2 && cd ../


