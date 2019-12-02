#pragma once
#include "zlib.h"
#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

#include "file_handling.hpp"
#include "load_fasta.hpp"
#include "nucl_conv.hpp"
#include "myargparser.hpp"
#include "constants.hpp"


using v_un_ch = std::vector<unsigned char>;


// bool check_file_exists(std::string filename);
// void filter_ref_sites(const std::string & selected_chrom, std::string & filename, char * ref);
// void filter_ref_bed(const std::string & selected_chrom, std::string & filename, char * ref);

/// STRUCTS
struct my_cov_rg {
  size_t nocpg=0, cpg=0;
  // my_cov_rg(size_t & rgs_size){
  //   nocpg.resize(rgs_size, 0);
  //   cpg.resize(rgs_size, 0);
  // }
  // std::vector<size_t> nocpg, cpg;
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


struct per_mle_run_backup {
  std::vector<size_t> idx_to_include; // , positions;
  size_t total_depth, min_pos, max_pos, curr_pos, n_cpgs;
  per_mle_run_backup (const pre_calc_per_site & d, const size_t & idx){
    n_cpgs = 1;
    idx_to_include.push_back(idx);
    // positions.push_back(d.position);
    total_depth = d.depth;
    min_pos = d.position;
    max_pos = d.position;
    curr_pos = d.position;
  }
} ;



struct alignment_data {
  alignment_data () {
    t_seq.reserve(100), t_posi.reserve(100), t_isop.reserve(100), t_qs.reserve(100), t_ref.reserve(100), t_positions.reserve(100);
  }
  unint strand, mapQ, n_nucleotides;
  std::vector<unint> t_seq, t_posi, t_isop, t_qs, t_ref, t_positions;

};


// NEW structs

using unint = unsigned int;

struct ObsF {
  ObsF(const double &_M,
       const double &_noM,
       const double &_me
       ):
    M(_M), noM(_noM), me(_me){}
  double M, noM;
  double me;
};

struct Obs {
  Obs(const unint &_pr, // prime 0,1
      const unint &_st, // strand 0,1
      const unint &_rp, // readpos int
      const unint &_rl,  // read length 0,1
      const unint &_bc, // base composition 0,1,2
      const unint &_se, // seq error 1-40
      const unint &_me // mapping error 1-40
      ):
    pr(_pr), st(_st),
    rp(_rp), rl(_rl),
    bc(_bc), se(_se),
    me(_me) {

  }

  unint pr : 1;  // 0..1  (1 bits)
  unint st : 1;  // 0..1  (1 bits)
  unint rp : 6;  // 0..64 (6 bits) 8
  unint rl : 1;  // 0..1  (1 bits)
  unint bc : 2;  // 0,1,2 (2 bits) 3
  unint se : 6;  // 0-40  (6 bits)
  unint me : 6;  // 0-40  (6 bits)
};

using uni_ptr_obs = std::vector<std::unique_ptr<Obs>>;


struct Site_s {
  Site_s(){
    remaining_dinucl_genotypes.resize(6, 0);
  };

  void load_data(const double &_M, const double &_noM, const double &_me){
    data.push_back(ObsF(_M, _noM, _me));
    depth++;
  };
  size_t pos, depth=0;
  std::vector<double> remaining_dinucl_genotypes;
  std::vector<ObsF> data;

};

struct per_mle_run {
  std::vector<std::unique_ptr<Site_s>> to_include; // , positions;
  std::vector<size_t> positions;
  Site_s center_site;
  size_t total_depth, min_pos, max_pos, curr_pos, n_cpgs;
  per_mle_run (const Site_s & _d){
    n_cpgs = 1;
    to_include.emplace_back(std::make_unique<Site_s>(_d));
    center_site = _d;
    positions.push_back(_d.pos);
    total_depth = _d.depth;
    min_pos = _d.pos;
    max_pos = _d.pos;
    curr_pos = _d.pos;
  };

  void update_left(const Site_s & _d){
    min_pos = _d.pos;
    stats_update(_d);
  }

  void update_right(const Site_s & _d){
    max_pos = _d.pos;
    stats_update(_d);
  }

  void stats_update(const Site_s & _d){
    to_include.emplace_back(std::make_unique<Site_s>(_d));
    positions.push_back(_d.pos);
    total_depth += _d.depth;
    n_cpgs++;
  }

} ;


struct rgs_info {
  bool rg_split=false;
  std::vector<std::string> rgs;
  std::vector<size_t> cycles;
  size_t n = 0;
};

struct deamrates_void {
  deamrates_void(general_settings * _settings,
                 uni_ptr_obs &_cpg_data,
                 uni_ptr_obs &_nocpg_data){
    settings = _settings;
    cpg_data = std::move(_cpg_data);
    nocpg_data = std::move(_nocpg_data);
  };
  general_settings * settings;
  uni_ptr_obs cpg_data, nocpg_data;
  size_t iteration=0;

};


struct F_void {
  general_settings * settings;
  std::vector<std::unique_ptr<Site_s>> to_include;
  size_t iteration;
  F_void(general_settings * _s, std::vector<std::unique_ptr<Site_s>> & _m){
    settings = _s;
    to_include = std::move(_m);
    iteration = 0;
  }
};



// FUNCTIONS


inline double oplusnatl(const double & x, const double & y );

template <typename T>
inline T oplusInitnatl(const T & x,const T & y );

// template <class T>
// void checkfilehandle(T &fh, std::string filename);

// bool check_file_exists(std::string filename);

// double phred_to_double(int & phred);

// std::vector<double> phred_to_double_converter();


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
void get_pos_and_prime(general_settings & settings, const size_t & dist_3p, const size_t & dist_5p, size_t & pos, size_t & prime);
bool nucleotide_missing(const int & r_base1, const int & r_base2, const int & s_base1, const int & s_base2);
bool base_is_indel(const int & base1);
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
