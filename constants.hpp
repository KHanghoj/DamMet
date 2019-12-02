#pragma once

#include <vector>
#include <string>

const size_t PRIMES=2;
const size_t STRANDS=2;
const size_t NUCLEOTIDES=4;
const size_t METHSTATES=2;
const size_t METHSTATE=0;
const size_t UNMETHSTATE=1;

const size_t READPOS_MULT = PRIMES*STRANDS*std::pow(NUCLEOTIDES,4);
const size_t PRIME_MULT = STRANDS*std::pow(NUCLEOTIDES,4);
const size_t STRAND_MULT = std::pow(NUCLEOTIDES,4);
const size_t B1_MULT = std::pow(NUCLEOTIDES,3);
const size_t B2_MULT = std::pow(NUCLEOTIDES,2);
const size_t B3_MULT = NUCLEOTIDES;

const size_t N_types_TM = 7;

const size_t MAX_PHRED = 255;

const double SMALLTOLERANCE = 1e-7;

// FIXME: should reflect the presence of the observed base
double BASE_FREQ_FLAT_PRIOR = 0.25;
double DINUCL_FLAT_PRIOR = 1.0/7.0;

std::string UNKNOWN_RG("UNKNOWN");
std::string ALL_RG("ALL_DAMMET");
size_t ALL_DAMMET_RG_IDX = 0;

size_t do_haploid_model = 0;


using unint = unsigned int;
