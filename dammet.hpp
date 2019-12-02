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
#include "structs.hpp"

using v_un_ch = std::vector<unsigned char>;


// bool check_file_exists(std::string filename);
// void filter_ref_sites(const std::string & selected_chrom, std::string & filename, char * ref);
// void filter_ref_bed(const std::string & selected_chrom, std::string & filename, char * ref);



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
