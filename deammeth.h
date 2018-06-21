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
using dinucl_pair_of_pair = std::pair<
  std::pair<int, int>, // first genotype
  std::pair<int, int> // second genotype
  >;

using dinucl_pair_of_pairs = std::vector<dinucl_pair_of_pair>;


const double SMALLTOLERANCE = 1e-7;
const size_t METHSTATES=2;
const size_t METHSTATE=0;
const size_t UNMETHSTATE=1;
const size_t PRIMES=2;
const size_t STRANDS=2;
const size_t NUCLEOTIDES=4;

// tmf: // readpos, prime, strand, refdinucl1, refdinucl2, sampledinucl1,sampledinucl2
// FIXME: should reflect the presence of the observed base
double BASE_FREQ_FLAT_PRIOR = 0.25;

const std::vector<double> get_log_prior(){
  std::vector<double> res(7, std::log(0));
  res[0] = std::log(1);
  for (const auto & val : res ){
    std::cerr << val << '\n';
  }
  return res;
}

const std::vector<double> get_log_prior_type_specific(){
  std::cerr << "TYPE SPECIFIC" <<  std::endl;
  std::vector<double> res(7, 0);
  res[1] = 0.05;
  res[2] = res[1];
  res[3] = 0.001;
  res[4] = res[3];
  res[5] = res[3];
  res[6] = res[3];
  double s=0;
  for (auto & val : res){
    s+=val;
    val=std::log(val);
  }
  res[0] = std::log(1-s);
  return res;
}

const std::vector<double> get_log_prior_flat(){
  std::vector<double> res(7, std::log(0.001/6));
  res[0] = std::log(1-0.001);
  return res;
}

const std::vector<double> LOG_PRIORS=get_log_prior();
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


// as i provide rgs.size() this has to be a 'const type & val' or just 'type val'
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
  std::vector<double> remaining_dinucl_genotypes;

  pre_calc_per_site(const per_site & d){
    summary.resize(3);
    depth=d.depth;
    position=d.position;
    pre_noM.reserve(depth);
    pre_M.reserve(depth);
    maperrors.reserve(depth);
    remaining_dinucl_genotypes.reserve(6);
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
