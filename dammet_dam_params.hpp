#pragma once

#include <algorithm> // find
#include <cmath>
#include <fstream>  // open file
#include <iostream> // stdout/stdin/stderr
#include <random> // std::normal_distribution;
#include <signal.h>  // catching signal interrupt
#include <sstream> //stringstream
#include <string>
#include <unordered_map>
#include <utility> // pair
#include <vector>

#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

#include "load_fasta.hpp"
#include "nucl_conv.hpp"
#include "myargparser.hpp"

// FIXME: should reflect the presence of the observed base
double BASE_FREQ_FLAT_PRIOR = 0.25;
double DINUCL_FLAT_PRIOR = 1.0/7.0;



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



/// STRUCTS


struct my_cov_rg {
  std::vector<size_t> nocpg, cpg;
  my_cov_rg(size_t rgs_size){
    nocpg.resize(rgs_size,0);
    cpg.resize(rgs_size,0);
  }
} ;

struct per_site {
  size_t position, depth, no_deam_CT_GA;
  size_t dinucl_max_idx;
  bool exclude_for_deam;
  std::vector<int> prime, strand, pos_to_end;
  std::vector<std::pair<int, int>> bases;
  std::vector<std::pair<double, double>> quals;
  std::vector<size_t> readlengths;
  std::vector<size_t> base_compos;
  std::vector<double> seqerrors;
  std::vector<double> maperrors;
  std::vector<size_t> rgs;
  per_site(){
    depth=0;
    no_deam_CT_GA=0;
    exclude_for_deam=false;
  }
} ;

struct per_site_nocpg {
  size_t depth;
  std::vector<int> prime, pos_to_end;
  std::vector<size_t> base_compos;
  std::vector<double> seqerrors;
  std::vector<double> maperrors;
  std::vector<size_t> readlengths;
  std::vector<size_t> rgs;
  per_site_nocpg(size_t & max_data_points_to_include){
    depth=0;
    prime.reserve(max_data_points_to_include);
    pos_to_end.reserve(max_data_points_to_include);
    base_compos.reserve(max_data_points_to_include);
    seqerrors.reserve(max_data_points_to_include);
    maperrors.reserve(max_data_points_to_include);
    std::vector<size_t> readlengths;
    rgs.reserve(max_data_points_to_include);
  }
} ;

struct pre_calc_per_site {
  size_t position, depth;
  std::vector<double> maperrors, pre_noM, pre_M;
  std::vector<double> remaining_dinucl_genotypes;

  pre_calc_per_site(const per_site & d){
    depth=d.depth;
    position=d.position;
    pre_noM.reserve(depth);
    pre_M.reserve(depth);
    maperrors.reserve(depth);
    remaining_dinucl_genotypes.reserve(6);
  }
} ;


struct per_mle_run {
  std::vector<size_t> idx_to_include, positions;
  size_t total_depth, min_pos, max_pos, curr_pos, n_cpgs;
  per_mle_run (const pre_calc_per_site & d, const size_t & idx){
    n_cpgs = 1;
    idx_to_include.push_back(idx);
    positions.push_back(d.position);
    total_depth = d.depth;
    min_pos = d.position;
    max_pos = d.position;
    curr_pos = d.position;
  }
} ;


struct alignment_data {
  size_t strand, mapQ, n_nucleotides;
  std::vector<size_t> t_seq, t_posi, t_isop, t_qs, t_ref, t_positions;
  alignment_data () {
    t_seq.reserve(100), t_posi.reserve(100), t_isop.reserve(100), t_qs.reserve(100), t_ref.reserve(100), t_positions.reserve(100);
  }
};



// FUNCTIONS

inline double oplusnatl(const double & x, const double & y );

template <typename T>
inline T oplusInitnatl(const T & x,const T & y );

template <typename T>
void checkfilehandle(T &fh, std::string filename);

bool check_file_exists(std::string filename);

double phred_to_double(int & phred);

std::vector<double> phred_to_double_converter();


//  generate long array with T template
template <typename T>
std::vector<T> init_tallymat(const size_t & readpos);
size_t get_idx_tm(const size_t & readpos,
                  const size_t & prime,
                  const size_t & strand,
                  const size_t & b1,
                  const size_t & b2,
                  const size_t & b3,
                  const size_t & b4);

// cpg positions per chrom
using keeplist_map = std::unordered_map<size_t, size_t>;
const keeplist_map get_cpg_genomic_pos(const char * ref, const size_t & seq_len);
// const std::vector<int> get_cpg_genomic_pos (const char * ref, const size_t & seq_len);
const std::vector<int> get_cpg_chrom_bool(const char * ref, const size_t & seq_le);



// from deammeth.h
using keeplist_map = std::unordered_map<size_t, size_t>;
// using keeplist_map = std::map<size_t, size_t>;


// header management:
// https://stackoverflow.com/a/2596554/2788987
using dinucl_pair_of_pair = std::pair<
  std::pair<int, int>, // first genotype
  std::pair<int, int> // second genotype
  >;

using dinucl_pair_of_pairs = std::vector<dinucl_pair_of_pair>;

std::vector<double> get_log_prior(){
  std::vector<double> res(7, std::log(0));
  res[0] = std::log(1);
  /* for (const auto & val : res ){ */
  /*   std::cerr << val << '\n'; */
  /* } */
  return res;
}

std::vector<double> get_log_prior_flat(){
  std::vector<double> res(7, std::log(0.001/6.0));
  res[0] = std::log(1-0.001);
  return res;
}

/* std::vector<double> get_log_prior_type_specific(){ */
/*   std::cerr << "TYPE SPECIFIC" <<  std::endl; */
/*   std::vector<double> res(7, 0); */
/*   res[1] = 0.05; */
/*   res[2] = res[1]; */
/*   res[3] = 0.001; */
/*   res[4] = res[3]; */
/*   res[5] = res[3]; */
/*   res[6] = res[3]; */
/*   double s=0; */
/*   for (auto & val : res){ */
/*     s+=val; */
/*     val=std::log(val); */
/*   } */
/*   res[0] = std::log(1-s); */
/*   return res; */
/* } */

std::vector<double> LOG_PRIORS=get_log_prior_flat();
// const std::vector<double> LOG_PRIORS=get_log_prior_flat();
// const std::vector<double> LOG_PRIORS=get_log_prior_type_specific();
// const std::vector<double> PRIORS=get_prior();



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
const dinucl_pair_of_pairs SEVEN_DINUCL_GENOTYPES = GENERATE_SEVEN_DINUCL_GENOTYPES();


// IO funcs
void dump_count_file(general_settings & settings, std::vector<int> & tm, std::string & rg);

std::vector<double> read_count_file(general_settings & settings, std::string & rg);


// alignment
void get_bases(const alignment_data & d, const int & b1, const int & b2, int & r_base1, int & r_base2, int & s_base1, int & s_base2);
void get_pos_and_prime(general_settings & settings, const size_t & dist_3p, const size_t & dist_5p, size_t & pos, size_t & prime);
bool nucleotide_missing(const int & r_base1, const int & r_base2, const int & s_base1, const int & s_base2);
bool base_is_indel(const int & base1, const int & base2);
alignment_data align_read(bam1_t * rd, char * ref);
bool no_ref_cpg(const int & r_base1, const int & r_base2);
bool cpg(const int & base1, const int & base2);

void add_aligned_data(general_settings &settings,
                      const alignment_data &d,
                      const keeplist_map &cpg_map,
                      const char *ref,
                      std::vector<per_site> &data,
                      std::vector<int> &tm,
                      per_site_nocpg &nocpg_data,
                      my_cov_rg & cov_rg,
                      size_t & rg_idx,
                      size_t & ncycles,
                      size_t & exclude_CnonCpG);




// as i provide rgs.size() this has to be a 'const type & val' or just 'type val'
