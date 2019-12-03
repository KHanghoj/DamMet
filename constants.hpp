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

// for random sampling with replacement
std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 rn_generator(rd()); //Standard mersenne_twister_engine seeded with rd()



using unint = unsigned int;

// header management:
// https://stackoverflow.com/a/2596554/2788987
using dinucl_pair_of_pair = std::pair<
  std::pair<int, int>, // first genotype
  std::pair<int, int> // second genotype
  >;
using dinucl_pair_of_pairs = std::vector<dinucl_pair_of_pair>;

const dinucl_pair_of_pairs GENERATE_SEVEN_DINUCL_GENOTYPES(){
  dinucl_pair_of_pairs res;

  for (int i=0; i<7; i++){
    dinucl_pair_of_pair di;
    di.first.first=1;
    di.first.second=1;
    di.second.first=2;
    di.second.second=2;
    res.push_back(di);
  }
  // see latex for explanation
  // related to methylation
  res[1].first.second=3;
  res[2].second.second=0;
  // unrelated to methylation:
  res[3].first.second=0;
  res[4].second.second=1;
  res[5].first.second=2;
  res[6].second.second=3;
  return res;
}


std::vector<double> GET_LOG_PRIOR(){
  std::vector<double> res(7, std::log(0));
  res[0] = std::log(1);
  return res;
}

std::vector<double> GET_LOG_PRIOR_FLAT(){
  std::vector<double> res(7, std::log(0.001/6.0));
  res[0] = std::log(1-0.001);
  return res;
}


double phred_to_double(size_t & phred){
  return std::pow(10, -((double)phred/10.0));
}
std::vector<double> phred_to_double_converter(){
  std::vector<double> res;
  for (size_t phred=0; phred<=MAX_PHRED; phred++){
    res.push_back(phred_to_double(phred));
  }
  return res;
}
