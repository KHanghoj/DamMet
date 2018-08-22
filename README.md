
# Deammeth, a full probabilistic model for mapping ancient methylomes
=====================================================================
Deammeth is probabilistic model for mapping ancient methylomes using sequencing data underlying an ancient specimen.
The model is implemented as a two step procedure. The first step recovers a maximum likelihood estimate (MLE) of the position specific deamination rates for methylated and unmethylated cytosine residues. The second step, making use of these deamination rates, returns a MLE of the methylation level in a user-defined genomic window. The two step procedure as implemented in deammeth is fully 

## Requirements
deammeth requires [NLopt](https://nlopt.readthedocs.io/en/latest/) to be installed.

    apt-get install nlopt-devel

Additionally, it makes use of [htslib](https://github.com/samtools/htslib.git) to parse genome references and BAM files. htslib requires [zlib](https://zlib.net/) to be installed.
## Installation

``` bash
  git clone https://gitlab.com/KHanghoj/deammeth.git
  cd deammeth
  git clone https://github.com/samtools/htslib.git; cd htslib ; make ; cd ..
  make
  cd ..
```

## Tutorial

## Simulating sequence data with methylation specific deamination patterns using [Gargammel](https://github.com/grenaud/gargammel).

## Citation

asdf asdf 
git clone https://github.com/samtools/htslib.git

g++  -O3  deammeth.cc myargparser.cc ./htslib/libhts.a -std=c++11 -I./htslib/ -lz -lm -lbz2 -llzma -lnlopt -o deammeth -lpthread
