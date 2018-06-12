DINGDONG

git clone https://github.com/samtools/htslib.git

g++  -O3  deammeth.cc myargparser.cc ./htslib/libhts.a -std=c++11 -I./htslib/ -lz -lm -lbz2 -llzma -lnlopt_cxx -o deammeth -lpthread
