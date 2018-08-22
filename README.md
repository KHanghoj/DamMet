# Deammeth, a full probabilistic model for mapping ancient methylomes #

-------------------------------------------------------------------------------

Deammeth is probabilistic model for mapping ancient methylomes using sequencing data underlying an ancient specimen.
The model is implemented as a two step procedure. The first step recovers a maximum likelihood estimate (MLE) of the position specific deamination rates for methylated and unmethylated cytosine residues. The second step, making use of these deamination rates, returns a MLE of the methylation level in a user-defined genomic window. The two step procedure as implemented in deammeth is fully 

## Requirements ##

deammeth requires [NLopt](https://nlopt.readthedocs.io/en/latest/) to be installed.

``` bash
    apt-get install nlopt-devel
```

Additionally, it makes use of [htslib](https://github.com/samtools/htslib.git) to parse genome references and BAM files. htslib requires [zlib](https://zlib.net/) to be installed.

## Installation ##

``` bash
  git clone https://gitlab.com/KHanghoj/deammeth.git
  cd deammeth
  git clone https://github.com/samtools/htslib.git; cd htslib ; make -j2 ; cd ..
  make
  cd ..
```

## How to run deammeth ##

deammeth takes three required arguments, a bam file (-b), a reference genome (-r), and the chromosome of interest (-c). To demonstrate that deammeth works and produces the expected output, we have made a small running example.

### Running example ###

``` bash
  git clone https://gitlab.com/KHanghoj/deammeth-tutorial.git
  cd deammeth-tutorial
  git clone https://gitlab.com/KHanghoj/deammeth.git
  cd deammeth/
  git clone https://github.com/samtools/htslib.git; cd htslib ; make -j2 ; cd ..
  make
  cd ..
  bash run.sh
```

This example also shows the several of the available options. The main script in this example ('run.sh') can easily be modified for analyses of any sample of interest.


### Options and input formats ###

All avaliable options with a description can displayed by running deammeth without any arguments (`./deammeth/deammeth`). 

#### Input formats ####

1. *-R* allows the user to provide a list of read groups (ID) that should be considered individually for estimating deamination rates. Readgroups not present in the input file will be ignored and their sequencing data is not considered for any analyses. The input file should contain a single read group name (ID) per line. If the end-user does not provide a file to *-R* all, all sequencing data with be merged into a single read group. The latter is the default in deammeth.
2. *-E* allows the user to provide genomic sites that should be excluded. The file show contain a site per line (e.g. chr20 100001). The genomic position should be 1-based.
3. *-e* allows the user to provide genomic regions that should be excluded. The file should take the form of a standard BED file.

## Simulating sequence data with methylation specific deamination patterns using [gargammel](https://github.com/grenaud/gargammel). ##

Along with the publication of deammeth, we also developed a new feature to [gargammel]()(https://github.com/grenaud/gargammel) that enables the end-user to simulate ancient DNA sequences with methylation specific deamination patterns. With this new feature, it has become possible to answer questions like, what is the accuracy with the current sequecing effort, what is my false discovery rate, what is the optimal trade off between genomic resolution (e.g. genomic window) and sequencing effort. 

## Citation ##

WILL COME.



