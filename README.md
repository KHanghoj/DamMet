[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/dammet/README.html)

# DamMet, a full probabilistic model for mapping ancient methylomes #

-------------------------------------------------------------------------------

DamMet is probabilistic model for mapping ancient methylomes using sequencing data underlying an ancient specimen.
The model is implemented as a two step procedure. The first step recovers a maximum likelihood estimate (MLE) of the position specific deamination rates for methylated and unmethylated cytosine residues. The second step, making use of these deamination rates, returns a MLE of the methylation level in a user-defined genomic window. The software is multithreaded

## Installation ##
DamMet is dependent on [htslib](https://github.com/samtools/htslib.git) and [nlopt](https://nlopt.readthedocs.io/en/latest). Both software will be downloaded and installed with DamMet.

``` bash
  git clone --recurse-submodules https://github.com/KHanghoj/DamMet.git
  ## or update existing git: git submodule update --init --recursive 
  cd DamMet && make && cd ..
```
[zlib](https://zlib.net/) and [cmake](https://cmake.org/download) should be globally installed.

## How to run DamMet ##

DamMet requires a single positional argument (estDEAM or estF). estDEAM will estimate deamination profiles and estF estimates methylation levels. 

In addition DamMet also requires three keyword arguments:
- BAM file (-b)
- Reference genome (-r)
- a comma separated list of chromosomes/contigs (-c) or path to file with chromosomes/contigs (-cf).


The two step procedure could look like this:

```bash
    ./DamMet estDEAM -b ${BAM} -r ${FASTA} -c ${CHROM} -M ${M} -O ${OUT_PREFIX}
    ./DamMet estF -b ${BAM} -r ${FASTA} -c ${CHROM} -M ${M} -O ${OUT_PREFIX} -N ${NCPG}
```

To obtain read depth for each CpG use the following command:
```bash
    ./DamMet getSites -b ${BAM} -r ${FASTA} -c ${CHROM} -O ${OUT_PREFIX}
```
This can be useful to ensure that the exact same CpGs are considered in a cross sample comparison.


### Running example ###

To demonstrate that DamMet works and produces the expected output, we have made a small running example.

``` bash
git clone https://github.com/KHanghoj/DamMet.git
make -C DamMet 
git clone https://github.com/KHanghoj/DamMet-tutorial.git
cd DamMet-tutorial
bash run.sh
```

The output plot (*result.pdf*) of this running example should be identical to *result.expected.pdf*. The main script in this example ('run.sh') can easily be used as a template to suit any analyses of interest.


### Options, Input and Output formats ###

All available options followed by a description can displayed by running DamMet without any arguments (`./DamMet/DamMet`). 

#### Special Input formats ####

1. *-R* allows the user to provide a list of read groups (ID) that should be considered individually for estimating deamination rates. Readgroups not present in the input file will be ignored and their sequencing data is not considered for any analyses. The input file should contain a single read group name (ID) per line. If the user does not provide a file to *-R* all, all sequencing data with be merged into a single read group. The latter is the default in DamMet.
2. *-E* allows the user to provide genomic sites that should be excluded. The file show contain a site per line (e.g. chr20 100001). The genomic position should be 1-based.
3. *-e* allows the user to provide genomic regions that should be excluded. The file should take the form of a standard BED file.

#### Output format ####

Two types of files are produced by DamMet, namely a *CHR.READGROUP.deamrates* file and a *CHR.READGROUP.[BED].F* file.

##### *CHR.READGROUP.deamrates* #####
This file contains the MLE of the deamination rates. 
1. Methylation status (Methylated == 0; UnMethylated == 1)
2. DNA position along a DNA molecule (0-based)
3. Prime (5-prime == 0; 3-prime == 1)
4. Deamination rate

##### *CHR.READGROUP.CHROM.BED.F* or *CHR.READGROUP.CHROM.F* #####
Every row is a BED region or a genomic CpG. Every column output comes with a description.


## Simulating sequence data with methylation specific deamination patterns using [gargammel](https://github.com/grenaud/gargammel). ##


Along with the publication of DamMet, we also developed a new feature to [gargammel](https://github.com/grenaud/gargammel) that enables the user to simulate ancient DNA sequences with methylation specific deamination patterns. With this new feature, it has become possible to answer questions like, what is the accuracy with the current sequencing effort and what is the optimal trade off between genomic resolution (e.g. genomic window) and sequencing effort. 

## Citation ##

https://academic.oup.com/gigascience/article/8/4/giz025/5475519

# Troubleshooting #
## htslib installation issues ##
For one user installing the following software fixed some issues installing htslib:
`libbz2-1.0 libbz2-dev libbz2-ocaml libbz2-ocaml-dev liblzma-dev`.
