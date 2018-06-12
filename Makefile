deammeth: deammeth.cc myargparser.cc
	g++  -O3  deammeth.cc myargparser.cc ./htslib/libhts.a -std=c++11 -I./htslib/ -lz -lm -lbz2 -llzma -lnlopt_cxx -o deammeth -lpthread
