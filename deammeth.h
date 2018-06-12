#include <signal.h>  // catching signal interrupt
#include <algorithm> // find
#include <cmath>
#include <utility> // pair
#include <unordered_map>
#include <fstream>
#include <sstream> //stringstream
#include "load_fasta.h"
#include "myargparser.h"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <random> // std::normal_distribution;

using keeplist_map = std::unordered_map<size_t, size_t>;
// using keeplist_map = std::map<size_t, size_t>;


// header management:
// https://stackoverflow.com/a/2596554/2788987
using dinucl_pair_of_pair = std::pair<std::pair<int, int>, std::pair<int, int> >;
using dinucl_pair_of_pairs = std::vector<dinucl_pair_of_pair>;
using dinucl_pair = std::pair<int, int>;
using dinucl_pairs = std::vector<dinucl_pair>;

const double SMALLTOLERANCE = 1e-7;
const size_t METHSTATES=2;
const size_t METHSTATE=0;
const size_t UNMETHSTATE=1;
const size_t PRIMES=2;
const size_t STRANDS=2;
const size_t NUCLEOTIDES=4;

// tmf: // readpos, prime, strand, refdinucl1, refdinucl2, sampledinucl1,sampledinucl2

// FIXME: should be user specificed options
double MIN_DINUCL_GENOTYPE_PROB = 0.0001;

// FIXME: should reflect the presence of the observed base
double BASE_FREQ_FLAT_PRIOR = 0.25;

const dinucl_pairs GENERATE_DINUCL_PAIRS(){
  dinucl_pairs res;
  for (int i=0; i<4; i++){
      for (int j=0; j<4; j++){
        dinucl_pair temp (i,j);
        res.push_back(temp);
      }
  }
  return res;
}

const dinucl_pairs ALL_DINUCL_PAIRS = GENERATE_DINUCL_PAIRS();

dinucl_pair_of_pairs GENERATE_DINUCL_PAIR_OF_PAIRS(){
  dinucl_pair_of_pairs res;
  res.reserve(ALL_DINUCL_PAIRS.size()*ALL_DINUCL_PAIRS.size());
  for (auto it1=ALL_DINUCL_PAIRS.begin(); it1!=ALL_DINUCL_PAIRS.end(); it1++){
    for (auto it2=it1; it2!=ALL_DINUCL_PAIRS.end(); it2++){
      dinucl_pair_of_pair temp_res;
      temp_res.first = *it1;
      temp_res.second = *it2;
      res.push_back(temp_res);
    }
  }
  return res;
}

size_t find_pair_pair_idx(const int & r_base1, const int & r_base2, const int & s_base1, const int & s_base2){
  size_t counter =0;
  bool found = false;
  for (const auto & pairs : GENERATE_DINUCL_PAIR_OF_PAIRS()) {
    if (pairs.first.first == r_base1 && pairs.first.second == r_base2 && pairs.second.first == s_base1 && pairs.second.second == s_base2 ){
      found = true;
      break;
    }
    counter++;
  }
  if(!found){
    std::cerr << "Could not find IDX " << r_base1 << r_base2 << " " << s_base1 << s_base2 << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
  return counter;
}

dinucl_pair_of_pairs ALL_DINUCL_PAIR_OF_PAIRS =  GENERATE_DINUCL_PAIR_OF_PAIRS();


size_t IDX_CpG_CpG = find_pair_pair_idx(1,2,1,2);  // 81
size_t IDX_CpG_CpA = find_pair_pair_idx(1,0,1,2);  // 60
size_t IDX_CpG_TpG = find_pair_pair_idx(1,2,3,2);  // 89
size_t IDX_TpG_TpG = find_pair_pair_idx(3,2,3,2);  // 133
size_t IDX_CpA_CpA = find_pair_pair_idx(1,0,1,0);  // 58

// as i provide rgs.size() this has to be a 'const type & val' or just 'type val'
struct my_cov_rg {
  std::vector<size_t> nocpg, cpg;
  my_cov_rg(size_t rgs_size){
    nocpg.resize(rgs_size,0);
    cpg.resize(rgs_size,0);
  }
} ;

struct per_site {
  size_t position, depth;
  size_t dinucl_max_idx;
  std::vector<int> prime, strand, pos_to_end;
  std::vector<std::pair<int, int>> bases;
  std::vector<std::pair<int, int>> quals;
  std::vector<size_t> base_compos;
  std::vector<double> seqerrors;
  std::vector<double> maperrors;
  std::vector<size_t> rgs;
  per_site(){
    depth=0;
  }
} ;

struct per_site_nocpg {
  size_t depth;
  std::vector<int> prime, pos_to_end;
  std::vector<size_t> base_compos;
  std::vector<double> seqerrors;
  std::vector<double> maperrors;
  std::vector<size_t> rgs;
  per_site_nocpg(size_t & max_data_points_to_include){
    depth=0;
    prime.reserve(max_data_points_to_include);
    pos_to_end.reserve(max_data_points_to_include);
    base_compos.reserve(max_data_points_to_include);
    seqerrors.reserve(max_data_points_to_include);
    maperrors.reserve(max_data_points_to_include);
    rgs.reserve(max_data_points_to_include);
  }
} ;

struct pre_calc_per_site {
  size_t position, depth;
  std::vector<double> maperrors, pre_noM, pre_M, summary;
  pre_calc_per_site(const per_site & d){
    summary.resize(3);
    depth=d.depth;
    position=d.position;
    pre_noM.reserve(depth);
    pre_M.reserve(depth);
    maperrors.reserve(depth);
  }
} ;

struct pre_calc_per_site2 {
  size_t position, depth;
  std::vector<double> maperrors, pre_noM, pre_M, summary;
  pre_calc_per_site2(size_t & d, size_t & pos){
    summary.resize(3);
    depth=d;
    position=pos;
    pre_noM.reserve(depth);
    pre_M.reserve(depth);
    maperrors.reserve(depth);
  }
} ;

// struct per_em_run {
struct per_mle_run2 {
  std::vector<size_t> idx_to_include, summary, positions;
  size_t total_depth, min_pos, max_pos, curr_pos, n_cpgs;
  per_mle_run2 (const pre_calc_per_site2 & d, const size_t & idx){
    summary.resize(3);
    for(size_t i=0; i<d.summary.size(); i++){
      summary[i] = d.summary[i];
    }
    n_cpgs = 1;
    idx_to_include.push_back(idx);
    positions.push_back(d.position);
    total_depth = d.depth;
    min_pos = d.position;
    max_pos = d.position;
    curr_pos = d.position;
  }
} ;

struct per_mle_run {
  std::vector<size_t> idx_to_include, summary, positions;
  size_t total_depth, min_pos, max_pos, curr_pos, n_cpgs;
  per_mle_run (const pre_calc_per_site & d, const size_t & idx){
    summary.resize(3);
    for(size_t i=0; i<d.summary.size(); i++){
      summary[i] = d.summary[i];
    }
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
