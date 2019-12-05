#include <set>
#include <ctime>
#include <iomanip> // setprecision
#include <algorithm> // find
#include <cmath>
#include <fstream>  // open file
#include <iostream> // stdout/stdin/stderr
#include <random> // std::normal_distribution;
#include <sstream> //stringstream
#include <string>
#include <unordered_map>
#include <utility> // pair
#include <memory>
#include <vector>
#include <random>
#include <numeric> // accumulate
#include <thread>
#include <mutex>


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
#include "nlopt.hpp"

using v_un_ch = std::vector<unsigned char>;
using keeplist_map = std::unordered_map<size_t, size_t>;

bool VERBOSE = false;
int READS = 1e6;

std::mutex mtx;

std::vector<double> LOG_PRIORS=GET_LOG_PRIOR_FLAT();
// const std::vector<double> LOG_PRIORS=get_log_prior_flat();
// const std::vector<double> LOG_PRIORS=get_log_prior_type_specific();
// const std::vector<double> PRIORS=get_prior();
const dinucl_pair_of_pairs SEVEN_DINUCL_GENOTYPES = GENERATE_SEVEN_DINUCL_GENOTYPES();
std::vector<double> PHRED_TO_PROB_CONVERTER = phred_to_double_converter();


template <typename T>
void checkfilehandle(T &fh, std::string filename){
  if (! fh.is_open()){
    std::cerr << "Couldnt open file: " << filename << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
}


template <typename T>
void print_struct(std::unique_ptr<Obs> &data, T & out){
  out << data->pr << " " <<
    data->st << " " <<
    data->rp << " " <<
    data->rl << " " <<
    data->bc << " " <<
    data->se << " " <<
    data->me << '\n';
}

void print_data(uni_ptr_obs &cpg_data,
                uni_ptr_obs &nocpg_data
                ){
  for (auto &x: cpg_data){
    print_struct(x, std::cout);
  }

  for (auto &x: nocpg_data){
    print_struct(x, std::cout);
  }
}

void print_log(general_settings & settings, bool to_screen=true){
  if(to_screen)
    std::cerr << settings.buffer;
  settings.args_stream << settings.buffer;
  settings.buffer.clear();
}

// Returns logl( expl(x)+expl(y) )
inline double oplusnatl(const double & x, const double & y ){
    return x > y
         ? x + std::log1p(std::exp( y-x ) )
         : y + std::log1p(std::exp( x-y ) )  ;
}

// Returns log( expl(x)+expl(y) ), but does so without causing
// overflow or loss of precision.
// discards x if not initialized (x=0)
template <typename T>
inline T oplusInitnatl(const T & x,const T & y ){

    if( x == 0 ){ //no initialized, as log = 0 should not exist
        return y;
    }

    return x > y
      ? x + std::log1p( std::exp( y-x ) )
      : y + std::log1p( std::exp( x-y ) )  ;
}

const keeplist_map get_cpg_chrom_pos(const char * ref, const size_t & seq_len){
  keeplist_map res;
  size_t counter = 0;
  for (size_t i=0; i<seq_len-1; i++){
    if (refToInt[(int)ref[i]] == 1 && refToInt[(int)ref[i+1]] == 2){
      res[i] = counter;
      counter++;
    }
  }
  return res;
}

const v_un_ch get_cpg_chrom(const char * ref, const size_t & seq_len, int &ncpgs){
  v_un_ch res(seq_len, 0);
  for (size_t i=0; i<seq_len-1; i++){
    if ((refToInt[(int)ref[i]] == 1)&&(refToInt[(int)ref[i+1]] == 2)){
      res[i]=1;
      res[i+1]=1;
      ncpgs++;
    }
  }
  return res;
}


int my_divmod(const int & a, const int &b, int & res){
  res = a/b;
  return(a%b);
}

void idx_to_params_tm_fast(const int & idx, std::vector<int> & res){
  int rem;
  rem = my_divmod(idx, READPOS_MULT, res[0]);
  rem = my_divmod(rem, PRIME_MULT, res[1]);
  rem = my_divmod(rem, STRAND_MULT, res[2]);
  rem = my_divmod(rem, B1_MULT, res[3]);
  rem = my_divmod(rem, B2_MULT, res[4]);
  res[6] = my_divmod(rem, B3_MULT, res[5]);
}


template <typename T>
std::vector<T> init_tallymat(const size_t & readpos){
  int n = 1;
  // +1 to include the the array for alle obs within the reads
  n *= (readpos+2);
  n *= PRIMES * STRANDS;
  n *= std::pow(NUCLEOTIDES,4);
  return std::vector<T>(n, 0);
}

size_t get_idx_tm(const size_t & readpos,
                  const size_t & prime,
                  const size_t & strand,
                  const size_t & b1,
                  const size_t & b2,
                  const size_t & b3,
                  const size_t & b4){
  // https://en.wikipedia.org/wiki/Row-_and_column-major_order Address calculation in general
  return (readpos * READPOS_MULT +
          prime * PRIME_MULT +
          strand * STRAND_MULT +
          b1 * B1_MULT +
          b2 * B2_MULT +
          b3 * B3_MULT +
          b4);
}


void dump_count_file(general_settings & settings, std::vector<int> & tm, std::string & rg){
  std::string filename = settings.outbase + "." + rg + ".tallycounts";
  std::ofstream f (filename.c_str());
  checkfilehandle<std::ofstream>(f, filename);
  std::vector<int> info(N_types_TM, 0);
  if (f.is_open()) {
    for (int i=0; i<tm.size(); i++)
      if(tm[i]){
        idx_to_params_tm_fast(i, info);
        for (auto &x : info)
          f << x << ' ';
        f << tm[i] << '\n';
      }
  }
  f.flush();
  f.close();
}

std::vector<double> read_count_file(general_settings & settings, std::string & rg){
  std::vector<int> temp = init_tallymat<int>(settings.max_pos_to_end);
  std::vector<double> res = init_tallymat<double>(settings.max_pos_to_end);
  std::string filename = settings.outbase + "." + rg + ".tallycounts";
  std::ifstream f (filename.c_str());
  checkfilehandle<std::ifstream>(f, filename);

  if(f.is_open()){
    std::string row;
    size_t p, pr, st, r1, r2, s1, s2;
    double val;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> p >> pr >> st >> r1 >> r2 >> s1 >> s2 >> val;
      temp[get_idx_tm(p, pr, st, r1, r2, s1, s2)] = val;
      ss.clear();
    }
  }
  f.close();

  int idx;
  size_t readpos = settings.max_pos_to_end;
  for (size_t p = 0; p <= (readpos + 1); p++) {
    for (size_t pr = 0; pr < PRIMES; pr++) {
      for (size_t st = 0; st < STRANDS; st++) {
        for (size_t r1 = 0; r1 < NUCLEOTIDES; r1++) {
          for (size_t r2 = 0; r2 < NUCLEOTIDES; r2++) {
            size_t totalsum = 0;
            for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
              for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                idx = get_idx_tm(p, pr, st, r1, r2, s1, s2);
                if (temp[idx]) {
                  totalsum += temp[idx];
                }
              }
            }
            for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
              for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                idx = get_idx_tm(p, pr, st, r1, r2, s1, s2);
                if (temp[idx] && totalsum) {
                  res[idx] = (temp[idx] / (double)totalsum);
                }

                if (res[idx] < SMALLTOLERANCE){
                  res[idx] = SMALLTOLERANCE;
                }
                if(res[idx] > 1-SMALLTOLERANCE){
                  res[idx] = 1-SMALLTOLERANCE;
                }
              }
            }
          }
        }
      }
    }
  }
  return res;
};


void get_pos_and_prime(general_settings & settings, const size_t & dist_3p, const size_t & dist_5p, unint & pos, unint & prime){
    if ((dist_3p > settings.max_pos_to_end) &&
        (dist_5p > settings.max_pos_to_end)) {
      pos = settings.max_pos_to_end + 1;
      prime = 0; // 5p+1 by default
    } else if (dist_3p < dist_5p) {
      pos = dist_3p;
      prime = 1;
    } else {
      pos = dist_5p;
      prime = 0;
    }
}

bool base_is_indel(const int & b){
  return (b==5);
}

alignment_data align_read(bam1_t * rd, char * ref){

  alignment_data d ;
  d.n_nucleotides = rd->core.l_qseq;
  d.strand = bam_is_rev(rd);
  d.mapQ = rd->core.qual;
  if(d.mapQ>37){
    d.mapQ = 37;
  }
  // https://github.com/samtools/htslib/blob/bceff25af7bdf6a5dae77063d8719cdd92f56012/htslib/sam.h#L97
  uint8_t *quals = bam_get_qual(rd);
  uint8_t *seq = bam_get_seq(rd);
  int nCig = rd->core.n_cigar;
  uint32_t *cigs = bam_get_cigar(rd);
  int seq_pos = 0;              // position within sequence
  int wpos = rd->core.pos; // this value is the current position
                                // assocatied with the positions at seq_pos
  // see ANGSD code. pileup. heavily heavily inspired :)
  int hasInfo = 0;
  int c;
  int opCode, opLen;
  // loop through read by looping through the cigar string
  for (int i = 0; i < nCig; i++) {
    opCode = bam_cigar_op(cigs[i]);
    opLen = bam_cigar_oplen(cigs[i]);
    // opCode = cigs[i] & BAM_CIGAR_MASK;  // what to do
    // opLen = cigs[i] >> BAM_CIGAR_SHIFT; // length of what to do

    // fprintf(stderr,"opCode=%d opLen=%d seqPos=%d
    // wpos=%d\n",opCode,opLen,seq_pos,wpos);
    if (opCode == BAM_CINS ||
        opCode == BAM_CDEL) {             // handle insertions and deletions

      if (opCode == BAM_CINS && i == 0) { // skip indels if beginning of a read
        seq_pos += opLen;
        continue;
      }
      if (opCode == BAM_CINS && hasInfo == 0) {
        seq_pos += opLen;
        continue;
      }
      hasInfo++;
      if (opCode & BAM_CINS) {
        wpos--; // insertion/deletion is bound to the last position of the read
        for (int ii = 0; ii < opLen; ii++) {
          c = seq_nt16_int[bam_seqi(seq, seq_pos)];
          d.t_seq.push_back(c);
          // seq[ii] =  bam_is_rev(rd) ? tolower(c) : toupper(c);
          d.t_posi.push_back(seq_pos + 1);
          // posi[ii] = seq_pos + 1;
          d.t_isop.push_back(rd->core.l_qseq - seq_pos - 1);
          // isop[ii] = rd->core.l_qseq - seq_pos - 1;
          d.t_qs.push_back(quals[seq_pos]);
          // qs[ii] = quals[seq_pos];
          // 5 reflects an indel
          d.t_ref.push_back(5);

          d.t_positions.push_back(0);

          seq_pos++; // <- important, must be after macro
        }
        wpos++;
      } else { // this is the deletion part
        for(int ii=0;ii<opLen;ii++){
          // 5 reflects an indel
          d.t_seq.push_back(5);
          // seq[ii] =  bam_is_rev(rd) ? tolower(c) : toupper(c);
          d.t_posi.push_back(0);
          // posi[ii] = seq_pos + 1;
          d.t_isop.push_back(0);
          // isop[ii] = rd->core.l_qseq - seq_pos - 1;
          d.t_qs.push_back(0);
          // qs[ii] = quals[seq_pos];
          d.t_ref.push_back(refToInt[(int)ref[wpos +ii]]);

          d.t_positions.push_back(wpos+ii) ;

        }
        wpos += opLen;
      }

    }else if(opCode==BAM_CSOFT_CLIP){
      //occurs only at the ends of the read
      if(seq_pos == 0){
        //then we are at beginning of read and need to write mapQ
        seq_pos += opLen;
      }else//we are at the end of read, then break CIGAR loop
        break;

    }else if(opCode==BAM_CMATCH || opCode==BAM_CEQUAL || opCode==BAM_CDIFF) {

      hasInfo++;

      for(int fix=wpos ;wpos<(fix+opLen) ;wpos++) {
        // note that wpos is incrementing
        c = seq_nt16_int[bam_seqi(seq, seq_pos)];

        d.t_seq.push_back(c);
        // seq[ii] =  bam_is_rev(rd) ? tolower(c) : toupper(c);
        d.t_posi.push_back(seq_pos);
        // posi[ii] = seq_pos + 1;
        d.t_isop.push_back(rd->core.l_qseq - seq_pos - 1);
        // isop[ii] = rd->core.l_qseq - seq_pos - 1;
        d.t_qs.push_back(quals[seq_pos]);
        // qs[ii] = quals[seq_pos];
        d.t_ref.push_back(refToInt[(int)ref[wpos]]);

        d.t_positions.push_back(wpos) ;

        seq_pos++;
      }
    }else if(opCode==BAM_CREF_SKIP) {
      for(int ii=0;ii<opLen;ii++){
        d.t_seq.push_back(4);
        // seq[ii] =  bam_is_rev(rd) ? tolower(c) : toupper(c);
        d.t_posi.push_back(0);
        // posi[ii] = seq_pos + 1;
        d.t_isop.push_back(0);
        // isop[ii] = rd->core.l_qseq - seq_pos - 1;
        d.t_qs.push_back(0);
        // qs[ii] = quals[seq_pos];
        d.t_ref.push_back(refToInt[(int)ref[wpos +ii]]);
        d.t_positions.push_back(wpos+ii) ;
      }
      wpos += opLen;
      // dont care, i think
    }else if(opCode==BAM_CPAD||opCode==BAM_CHARD_CLIP) {
      //dont care
    }else{
      fprintf(stderr,"Problem with unsupported CIGAR opCode=%d\n",opCode);
    }
    //      exit(0);
  }
  //after end of read/parsing CIGAR always put the endline char


  if(d.t_seq.size() != d.t_posi.size()){
    std::cerr << "no good" << '\n';
    exit(EXIT_FAILURE);
  }

  return d;
}

bool no_ref_cpg(const int & r_base1, const int & r_base2){
  return(r_base1!=1 || r_base2!=2);
};

bool cpg(const int & base1, const int & base2){
  return(base1==1 && base2==2);
}

bool check_base_quality(const general_settings &settings,
                        const alignment_data &d,
                        const size_t &i_l, const size_t &i_r){
  bool trash=false;
  if (d.t_seq[i_l] >=4 || d.t_seq[i_r] >= 4 || // is read N
      d.t_ref[i_l] >=4 || d.t_ref[i_r] >= 4 || // is reference N
      d.t_qs[i_l] <settings.minbaseQ || d.t_qs[i_r] < settings.minbaseQ || // base quality
      (i_l>0 && (base_is_indel(d.t_ref[i_l-1]) || base_is_indel(d.t_seq[i_l-1]))) || // left is indel
      (i_r+1 <= d.t_seq.size()-1 && (base_is_indel(d.t_ref[i_r+1]) || base_is_indel(d.t_seq[i_r+1]))) // right is indel
      ) {
    trash = true;
  }
  return trash;
}

unint get_bc(const int & obs, // obs
             const int & no_ch, // ref
             const int & ch){ //  deamin/mut/error
  if(obs == no_ch){
    return 0;
  } else if (obs == ch){
    return 1;
  } else {
    return 2;
  }
}

size_t get_param_idx(const size_t & max_pos_to_end, const size_t & meth, const size_t & pos_to_end, const size_t & prime){
  return (meth * PRIMES * ( max_pos_to_end+2 ) + pos_to_end * PRIMES  + prime);
}

double base_condition_one_genotype(const unint &strand, const double &deam,
                                   const double &seqerror, const unint &base,
                                   const size_t &geno) {
  // - (strand==1 for reverse)
  // + (strand==0 for forward)
  // A->0; C->1; G->2; T->3
  double res;
  if ((strand==1 && geno == 2 && base == 2) ||
      (strand==0 && geno == 1 && base == 1)) {
    // - 2 2 or + 1 1
    res = ((1 - deam) * (1 - seqerror)) + (deam * seqerror / 3) +
          (seqerror * seqerror / 3);

  } else if ((strand==1 && geno == 2 && base == 0) ||
             (strand==0 && geno == 1 && base == 3)) {
    // - 2 0 or + 1 3
    res = (deam * (1 - seqerror)) + ((1 - deam) * seqerror / 3) +
          (seqerror * seqerror / 3);
    // std::cerr << strand << " " << geno << " " << base << " " << " " << deam << " " << res  << std::endl;
  } else {
    // the rest
    if (geno == base) {
      res = (1 - seqerror) + (seqerror * seqerror / 3);
    } else {
      res = (seqerror / 3) + (seqerror * seqerror / 3);
    }
  }
  return res;
}

double base_condition_both_gts(const unint &strand, const double &dam,
                               const double &s1, const unint &b1,
                               const double &s2, const unint &b2,
                               const dinucl_pair_of_pair &g){
  double prob_g1 =
    0.5 * base_condition_one_genotype(strand, dam, s1,
                                      b1, g.first.first) +
    0.5 * base_condition_one_genotype(strand, dam, s1,
                                      b1, g.first.second);
  double prob_g2 =
    0.5 * base_condition_one_genotype(strand, dam, s2,
                                      b2, g.second.first) +
    0.5 * base_condition_one_genotype(strand, dam, s2,
                                      b2, g.second.second);
  return prob_g1 * prob_g2;
}

void calc_M_noM(general_settings &settings,
                const alignment_data &d,
                const keeplist_map &cpg_map,
                size_t &ncycles,
                my_cov_rg &cov_rg,
                std::vector<double> deam,
                std::vector<Site_s> &cpg_data){

  unint dist_5p, dist_3p;
  unint pos, prime;
  unint base1, base2;
  unint rl = (d.n_nucleotides == ncycles);
  double seqerror1, seqerror2;
  double mapq = PHRED_TO_PROB_CONVERTER[d.mapQ];
  size_t cpg_idx;
  size_t boi; // base of interest

  double  prob_g1, prob_g2;
  double deamin_methylated, deamin_unmethylated;
  double noM, M;

  Site_s *s;

  bool anydata;

  for (size_t i = 0; i < d.t_seq.size()-1; i++){
    anydata = false;
    if (d.strand){ // true for negative strand
      boi = d.t_seq.size()-i-1; // closest to 5 prime
      if(check_base_quality(settings, d, boi-1, boi))
        continue;

      if(no_ref_cpg(d.t_ref[boi-1], d.t_ref[boi]))
        continue;

      dist_5p = d.t_isop[boi];
      dist_3p = d.t_posi[boi];
      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

      // if read == number of cycles. use the postions relative to the 5 prime.
      if (prime==1 && rl)
        get_pos_and_prime(settings, dist_5p*10, dist_5p, pos, prime);


      seqerror1 = PHRED_TO_PROB_CONVERTER[d.t_qs[boi-1]];
      seqerror2 = PHRED_TO_PROB_CONVERTER[d.t_qs[boi]];
      base1 = d.t_seq[boi-1];
      base2 = d.t_seq[boi];

#if 0
      std::cerr << d.t_positions[boi-1] << '\n';
#endif

      cpg_idx = cpg_map.at(d.t_positions[boi-1]);
      s = &cpg_data[cpg_idx];
      cov_rg.cpg++;
      anydata=true;
    } else {
      boi = i; // closest to 5 prime
      if(check_base_quality(settings, d, boi, boi+1))
        continue;

      if(no_ref_cpg(d.t_ref[boi], d.t_ref[boi+1]))
        continue;

      dist_5p = d.t_posi[boi];
      dist_3p = d.t_isop[boi];
      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

      // if read == number of cycles. use the postions relative to the 5 prime.
      if (prime==1 && rl)
        get_pos_and_prime(settings, dist_5p*10, dist_5p, pos, prime);

      seqerror1 = PHRED_TO_PROB_CONVERTER[d.t_qs[boi]];
      seqerror2 = PHRED_TO_PROB_CONVERTER[d.t_qs[boi+1]];
      base1 = d.t_seq[boi];
      base2 = d.t_seq[boi+1];

      cpg_idx = cpg_map.at(d.t_positions[boi]);
      s = &cpg_data[cpg_idx];
      cov_rg.cpg++;
      anydata=true;
    }

    if(anydata){
#if 0
      std::cerr << cpg_idx <<  " dingdong " << cov_rg.cpg << '\n';
#endif

#if 0
      std::cerr << pos << " " << prime << " " << d.strand << " " <<
        seqerror1 << " " << base1 << " " <<
        seqerror2 << " " << base2 << "\n";
#endif


      deamin_unmethylated = deam[get_param_idx(settings.max_pos_to_end, UNMETHSTATE, pos, prime)];
      noM = base_condition_both_gts(d.strand,
                                    deamin_unmethylated,
                                    seqerror1, base1,
                                    seqerror2, base2,
                                    SEVEN_DINUCL_GENOTYPES[0]);

      deamin_methylated = deam[get_param_idx(settings.max_pos_to_end, METHSTATE, pos, prime)];
      M = base_condition_both_gts(d.strand,
                                  deamin_methylated,
                                  seqerror1, base1,
                                  seqerror2, base2,
                                  SEVEN_DINUCL_GENOTYPES[0]);

      s->load_data(M, noM, mapq);
#if 0
      std::cerr << M << " " << noM << " " << s->depth << '\n';
#endif

      for (size_t dinucl_idx=1; dinucl_idx<SEVEN_DINUCL_GENOTYPES.size(); dinucl_idx++){
        double gt_it = base_condition_both_gts(d.strand,
                                               deamin_unmethylated,
                                               seqerror1, base1,
                                               seqerror2, base2,
                                               SEVEN_DINUCL_GENOTYPES[dinucl_idx]);

        s->remaining_dinucl_genotypes[dinucl_idx-1] +=
          std::log((1-mapq) * (gt_it) + mapq * DINUCL_FLAT_PRIOR);
#if 0
        std::cerr << s->remaining_dinucl_genotypes[dinucl_idx-1] << "\n";
#endif
      }
    }
  }
}

void add_aligned_data(general_settings &settings,
                      const alignment_data &d,
                      uni_ptr_obs &cpg_data,
                      uni_ptr_obs &nocpg_data,
                      std::vector<int> &tm,
                      my_cov_rg &cov_rg,
                      size_t &ncycles) {

  unint dist_5p, dist_3p;
  unint pos, prime, bc;
  unint rl = (d.n_nucleotides == ncycles);
  double seqerror = 0;
  size_t boi; // base of interest
  for (size_t i = 0; i < d.t_seq.size()-1; i++){

    if (d.strand){ // true for negative strand
      boi = d.t_seq.size()-i-1; // closest to 5 prime
      if(check_base_quality(settings, d, boi-1, boi))
        continue;

      if(d.t_ref[boi] != 2)
        continue;

      dist_5p = d.t_isop[boi];
      dist_3p = d.t_posi[boi];
      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

      // if read == number of cycles. use the postions relative to the 5 prime.
      if (prime==1 && rl)
        get_pos_and_prime(settings, dist_5p*10, dist_5p, pos, prime);
      bc = get_bc(d.t_seq[boi], 2, 0);

      if(no_ref_cpg(d.t_ref[boi-1], d.t_ref[boi]) &&
         cov_rg.nocpg < cov_rg.cpg){
        nocpg_data.emplace_back(std::make_unique<Obs>(prime,
                                                      d.strand,
                                                      pos,
                                                      rl,
                                                      bc,
                                                      d.t_qs[boi],
                                                      d.mapQ));
        cov_rg.nocpg++;
      }

      if(cpg(d.t_ref[boi-1], d.t_ref[boi])){
        cpg_data.emplace_back(std::make_unique<Obs>(prime,
                                                    d.strand,
                                                    pos,
                                                    rl,
                                                    bc,
                                                    d.t_qs[boi],
                                                    d.mapQ));
        cov_rg.cpg++;
      }
      tm[get_idx_tm(pos, prime, d.strand, d.t_ref[boi-1], d.t_ref[boi], d.t_seq[boi-1], d.t_seq[boi])]++;

    } else {
      boi = i; // closest to 5 prime
      if(check_base_quality(settings, d, boi, boi+1))
        continue;

      if(d.t_ref[boi] != 1)
        continue;

      dist_5p = d.t_posi[boi];
      dist_3p = d.t_isop[boi];
      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

      // if read == number of cycles. use the postions relative to the 5 prime.
      if (prime==1 && rl)
        get_pos_and_prime(settings, dist_5p*10, dist_5p, pos, prime);
      bc = get_bc(d.t_seq[boi], 1, 3);

      if(no_ref_cpg(d.t_ref[boi], d.t_ref[boi+1]) && cov_rg.nocpg < cov_rg.cpg){
          // adding data to nocpg
          nocpg_data.emplace_back(std::make_unique<Obs>(prime, d.strand, pos, rl, bc, d.t_qs[boi], d.mapQ));
          cov_rg.nocpg++;
      }

      if(cpg(d.t_ref[boi], d.t_ref[boi+1])){
        // adding data to cpg
        cpg_data.emplace_back(std::make_unique<Obs>(prime, d.strand, pos, rl, bc, d.t_qs[boi], d.mapQ));
        cov_rg.cpg++;
      }

      tm[get_idx_tm(pos, prime, d.strand, d.t_ref[boi], d.t_ref[boi+1], d.t_seq[boi], d.t_seq[boi+1])]++;
    }
  }
}


template <class T>
void print_single_array_parameters(const size_t & max_pos_to_end, std::vector<double> & param, T & out){
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==max_pos_to_end+1) || (pr==1 && p==0)){
          // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
        out << i << " " << p << " " << pr << " " << param[get_param_idx(max_pos_to_end, i, p, pr)] << '\n';
      }
    }
  }
}

double objective_func_deamrates(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
  if (!grad.empty()) {
    // resetting the array
    std::fill(grad.begin(), grad.end(), 0);
    // for (auto & val : grad){
    //   val = 0;
    // }
  }
  // https://stackoverflow.com/a/21249454/2788987

  deamrates_void *d = static_cast<deamrates_void*>(my_func_data);
  d->iteration++;

  double ll = 0;
  double res=0;
  double maperror, seqerror, deamin_methylated, deamin_unmethylated;
  double noM, M;
  double map_and_prior;
  size_t idx_param_meth, idx_param_unmeth;
  for (auto & p: d->cpg_data){

    seqerror = PHRED_TO_PROB_CONVERTER[p->se];
    maperror = PHRED_TO_PROB_CONVERTER[p->me];

    idx_param_meth = get_param_idx(d->settings->max_pos_to_end, METHSTATE, p->rp, p->pr);
    idx_param_unmeth = get_param_idx(d->settings->max_pos_to_end, UNMETHSTATE, p->rp, p->pr);
    deamin_methylated = x[idx_param_meth];
    deamin_unmethylated = x[idx_param_unmeth];
    map_and_prior =  maperror * BASE_FREQ_FLAT_PRIOR;
    res = 0;
    if (p->bc == 0) {
      // C->C -> (1-e) * (1-d) + (d * e/3)
      M = (1 - seqerror) * (1 - deamin_methylated) + (deamin_methylated * seqerror / 3.0);
      noM = (1 - seqerror) * (1 - deamin_unmethylated) + (deamin_unmethylated * seqerror / 3.0);
      res = (1-maperror) * (d->settings->M * M + (1-d->settings->M)*noM) + map_and_prior;

      if (!grad.empty()) {
        grad[idx_param_meth] += (d->settings->M*(1.0-maperror)*( 4.0*seqerror / 3.0 - 1) + map_and_prior) / res;
        grad[idx_param_unmeth] += ((1-d->settings->M) * (1.0-maperror)*( 4.0*seqerror / 3.0 - 1) + map_and_prior) /res;
      }

    } else if (p->bc == 1) {
      // C->T -> (1-e) * d + ((1-d) * e/3)
      M = (1 - seqerror) * (deamin_methylated) + ((1-deamin_methylated) * seqerror / 3.0);
      noM = (1 - seqerror) * deamin_unmethylated + ((1-deamin_unmethylated) * seqerror / 3.0);
      res = (1-maperror) * (d->settings->M * M + (1-d->settings->M)*noM)  + map_and_prior;
      if (!grad.empty()) {
        grad[idx_param_meth] += (d->settings->M * (1.0-maperror)*( -4.0*seqerror / 3.0 + 1) + map_and_prior) / res;
        grad[idx_param_unmeth] += ((1-d->settings->M) * (1.0-maperror)*( -4.0*seqerror / 3.0 + 1) + map_and_prior) / res;
      }
    } else {
      // C->[AG] -> (d)*e/3.0 + (1-d) * e/3.0
      M = (seqerror / 3) * deamin_methylated + ((1 - deamin_methylated) * seqerror / 3.0);
      noM = (seqerror / 3) * deamin_unmethylated + ((1 - deamin_unmethylated) * seqerror / 3.0);
      res = (1-maperror) * (d->settings->M * M + (1-d->settings->M)*noM)  + map_and_prior;

      // diff is zero
      // grad[idx_param_unmeth] += 0;
      // grad[idx_param_meth] += 0;

    }
    ll += std::log(res);
  }

  for (auto & p: d->nocpg_data){

    seqerror = PHRED_TO_PROB_CONVERTER[p->se];
    maperror = PHRED_TO_PROB_CONVERTER[p->me];
    idx_param_unmeth = get_param_idx(d->settings->max_pos_to_end, UNMETHSTATE, p->rp, p->pr);
    deamin_unmethylated = x[idx_param_unmeth];
    map_and_prior =  maperror * BASE_FREQ_FLAT_PRIOR;

   if (p->bc == 0) {
      // C->C -> (1-e) * (1-d) + (d * e/3)
      noM = (1 - seqerror) * (1 - deamin_unmethylated) + (deamin_unmethylated * seqerror / 3.0);
      res = (1-maperror) * noM  + map_and_prior;

      if (!grad.empty()) {
        grad[idx_param_unmeth] += ((1-maperror) * ( 4.0*seqerror / 3.0 - 1) + map_and_prior) / res;
      }

    } else if (p->bc == 1) {
      // C->T -> (1-e) * d + ((1-d) * e/3)
      noM = (1-seqerror) * deamin_unmethylated + ((1-deamin_unmethylated) * seqerror / 3.0);
      res = (1-maperror) * noM  + map_and_prior;
      if (!grad.empty()) {
        grad[idx_param_unmeth] += ((1-maperror) * ( -4.0*seqerror / 3.0 + 1) + map_and_prior )/ res;
      }
    } else {
      // C->[AG] -> (d)*e/3.0 + (1-d) * e/3.0
      noM = (seqerror / 3) * deamin_unmethylated + ((1 - deamin_unmethylated) * seqerror / 3.0);
      res = (1-maperror) * noM  + map_and_prior;

      // always zero
      // grad[idx_param_unmeth] += 0;

    }

    ll += std::log(res);
  }

  if(d->iteration%50==0){
    std::cerr << "\t-> Iteration: " << d->iteration << " llh: " << ll << '\r';
  }

  // if (!grad.empty()) {
  //   print_single_array_parameters(p->settings->max_pos_to_end, grad, std::cerr);
  // }
  // return -ll;
  return ll;
}

std::vector<double> single_array_parameters(const size_t & max_pos_to_end, const std::vector<double> & tmf){
  size_t max_length = PRIMES*(max_pos_to_end+2)*PRIMES ;
  std::vector<double> res(max_length, SMALLTOLERANCE);
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==max_pos_to_end+1) || (pr==1 && p==0)){
          // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }

        if(i==METHSTATE){ // methylated cpg
          res[get_param_idx(max_pos_to_end, i, p, pr)] = tmf[get_idx_tm(p, pr, 0, 1, 2, 3, 2)];
        } else {

          res[get_param_idx(max_pos_to_end, i, p, pr)] = tmf[get_idx_tm(p, pr, 0, 1, 1, 3, 1)];
        }
      }
    }
  }
  return res;
}

// -109750 with BOBYQA 10 minutes
// -109747 with MMA (requires derivatives) 20 seconds wuhu
// void run_deamrates_optim(general_settings & settings, const std::vector<double> & tmf, std::vector<per_site> data){

/// got this far
void run_deamrates_optim(general_settings & settings,
                         const std::vector<double> & tmf,
                         uni_ptr_obs &cpg_data,
                         uni_ptr_obs &nocpg_data,
                         std::string & rg){
  time_t start_time_optim,end_time_optim;
  time(&start_time_optim);
  std::string filename = settings.outbase + "." + rg + ".deamrates";
  std::ofstream f (filename.c_str());
  checkfilehandle<std::ofstream>(f, filename);
  std::vector<double> param = single_array_parameters(settings.max_pos_to_end, tmf);
  // print_single_array_parameters(settings.max_pos_to_end, param, std::cerr);

  deamrates_void void_stuff(&settings, cpg_data, nocpg_data);
  // print_single_array_parameters(settings.max_pos_to_end, param, std::cerr);
  std::vector<double> t; // just a dummy
  f << "##Initial (param are freq from tallycounts file) llh: " << objective_func_deamrates(param, t, &void_stuff) << std::endl;
  size_t n_params = PRIMES * ( settings.max_pos_to_end+2) * PRIMES ;
  // after 2000+ iterations (many minutes) : BOBYQA -239892
  // nlopt::opt opt(nlopt::LN_BOBYQA, n_params);
  // 100 iterations (60 seconds) : -239892
  nlopt::opt opt(nlopt::LD_MMA, n_params);

  // std::vector<double> lb(n_params, SMALLTOLERANCE), up(n_params, 1-SMALLTOLERANCE);
  std::vector<double> lb(n_params, SMALLTOLERANCE), up(n_params, 1-SMALLTOLERANCE);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(up);

  //opt.set_min_objective(objective_func_deamrates, &void_stuff);
  opt.set_max_objective(objective_func_deamrates, &void_stuff);

  opt.set_xtol_abs(0);
  opt.set_ftol_abs(1e-15);
  opt.set_xtol_rel(0);
  opt.set_ftol_rel(0);

  double minf;
  nlopt::result result = opt.optimize(param, minf);
  time(&end_time_optim);
  settings.buffer += "\t-> Iteration: " + std::to_string(void_stuff.iteration) +
    " llh: " + std::to_string(minf) + " in " +
    std::to_string(difftime(end_time_optim, start_time_optim)) + " seconds" + '\n';
  print_log(settings);

  f << "##Optim return code: " << result << std::endl;
  f << "##Iterations: " << void_stuff.iteration << std::endl;
  f << "##llh: " << minf << std::endl;

  print_single_array_parameters(settings.max_pos_to_end, param, f);
  //   print_single_array_parameters(settings.max_pos_to_end, param, std::cout);
  f.close();
}


std::vector<double> load_deamrates(general_settings & settings, std::string & rg){
  size_t max_length = PRIMES*(settings.max_pos_to_end+2)*PRIMES ;
  std::vector<double> res(max_length, SMALLTOLERANCE);
  std::string filename = settings.outbase + "." + rg + ".deamrates";
  std::ifstream f (filename.c_str());
  checkfilehandle<std::ifstream>(f, filename);
  if(f.is_open()){
    std::string row;
    size_t methstate, pos, prime;
    double rate;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      if(row.substr(0,1)=="#"){
        ss.clear();
        continue;
      }
      ss >> methstate >> pos >> prime >> rate;
      res[get_param_idx(settings.max_pos_to_end, methstate, pos, prime)] = rate;
      ss.clear();
    }
  }
  f.close();
  // print_single_array_parameters(settings.max_pos_to_end, res, std::cerr);
  return res;
}

std::vector<double> load_deamrates_f(general_settings & settings){
  size_t max_length = PRIMES*(settings.max_pos_to_end+2)*PRIMES ;
  std::vector<double> res(max_length, SMALLTOLERANCE);
  std::string filename = settings.deamrates_filename;
  std::ifstream f (filename.c_str());
  checkfilehandle<std::ifstream>(f, filename);
  if(f.is_open()){
    std::string row;
    size_t methstate, pos, prime;
    double rate;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      if(row.substr(0,1)=="#"){
        ss.clear();
        continue;
      }
      ss >> methstate >> pos >> prime >> rate;
      res[get_param_idx(settings.max_pos_to_end, methstate, pos, prime)] = rate;
      ss.clear();
    }
  }
  f.close();
  // print_single_array_parameters(settings.max_pos_to_end, res, std::cerr);
  return res;
}

// void update_mle_data(per_mle_run & mle_data, const pre_calc_per_site & data, size_t & idx){
//   mle_data.idx_to_include.push_back(idx);
//   mle_data.total_depth += data.depth;
//   mle_data.positions.push_back(data.position);
//   mle_data.n_cpgs++;
//   if (mle_data.min_pos > data.position) {
//     mle_data.min_pos = data.position;
//   } else if (mle_data.max_pos < data.position) {
//     mle_data.max_pos = data.position;
//   } else {
//     std::cerr << mle_data.curr_pos << " " << mle_data.min_pos << " " << mle_data.max_pos << " NOT GOOD POSITION " << data.position << '\n';
//   }
// }


double objective_func_F(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
  F_void *d = static_cast<F_void*>(my_func_data);
  d->iteration++;
  double log_like_cgcg_geno, read_der;
  double reads_grad, geno_grad;
  double read_like, log_like_genos ;
  double log_like_window=0;

  if(!grad.empty()){
    grad[0] = 0;
  }
  for (auto & site : d->to_include){
    log_like_cgcg_geno = 0;
    reads_grad = 0;
    geno_grad = 0;
    log_like_genos = 0;
    for (auto &obs : site->data){
      read_der = (1-obs.me) * (obs.M - obs.noM);
      read_like = (1-obs.me) * ((1-x[0]) * obs.noM + x[0] * obs.M) + obs.me*DINUCL_FLAT_PRIOR;
      log_like_cgcg_geno += std::log(read_like);
      if(!grad.empty()){
        // d/df per read
        reads_grad += read_der/read_like;
      }
    // wolfram: derivative  ln((1-w)*((1-x) * k + (x * h)) + w*p)
    // wolfram got it right, just modified a few things for readability
    // if(!grad.empty()){
    //   grad[0] += ((1-mape) * (M - noM)) / res;
    // }
    }

    if(!grad.empty()){
      // exp(log_like_cgcg_geno + LOG_PRIORS) * d/df all_reads. this is only including cgcg as the derivatives of the remaining are 0 as they are unrelated to F.
      geno_grad = std::exp(log_like_cgcg_geno + LOG_PRIORS[0]) * reads_grad;
    }

    log_like_genos = log_like_cgcg_geno + LOG_PRIORS[0];
    // summation of the likelihood of the remaining genotypes
    for (auto it=site->remaining_dinucl_genotypes.begin(); it!=site->remaining_dinucl_genotypes.end(); it++){
       log_like_genos = oplusnatl(log_like_genos, *it);
    }

    if(!grad.empty()){
      // convert log to raw with exp: 1/exp(ll) * d/df all_genos
      grad[0] +=  geno_grad / std::exp(log_like_genos);
    }
    log_like_window += log_like_genos;
  }
  return log_like_window;
}

double objective_func_F_second_deriv(const double & f,  void *my_func_data){
  F_void *d = static_cast<F_void*>(my_func_data);
  double total_d2=0;
  // double nominator, denominator;
  double noM, M, mape ;
  double log_like_cgcg_geno;
  double reads_grad, geno_grad;
  double read_der, read_like, log_like_genos ;
  double numerator_read_der2, denominator_read_der2, reads_der2;
  double geno_grad_der2, like_genos, like_genos_der2;
  for (auto & site : d->to_include){
    log_like_cgcg_geno = 0;
    reads_grad=0;
    reads_der2=0;
    for (auto &obs : site->data){
      read_der = (1-obs.me) * (obs.M - obs.noM);
      read_like = (1-obs.me) * ((1-f) * obs.noM + f * obs.M) + obs.me*DINUCL_FLAT_PRIOR;
      reads_grad += read_der/read_like;

      numerator_read_der2 = std::pow((1-obs.me),2) * std::pow((obs.M - obs.noM),2);
      denominator_read_der2 = std::pow((1-obs.me) * (obs.noM * (1-f) + f * obs.M) + obs.me*DINUCL_FLAT_PRIOR, 2);
      reads_der2 += - (numerator_read_der2 / denominator_read_der2);
      // \frac{\partial }{\partial \:x^2}\left(ln\left(\left(1-w\right)\cdot \left(\left(1-x\right)\:\cdot \:n\:+\:\left(x\:\cdot \:M\right)\right)\:+\:w\cdot \:p\right)\:\:\right)

      log_like_cgcg_geno += std::log(read_like);
      // https://www.symbolab.com/
      // wolfram: second derivative ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // second derivative:
      // second_der += (std::pow((1-obs.me),2) * std::pow((obs.M - obs.noM),2)) / std::pow((-obs.M * (obs.me-1) * f + obs.noM * (obs.me-1) * (1-f) + obs.me*DINUCL_FLAT_PRIOR), 2);
    }
    // exp(log_like_cgcg_geno + LOG_PRIORS) * d/df all_reads. this is only including cgcg as the derivatives of the remaining are 0 as they are unrelated to F.

    //g(x) first derivative     z(x)                              y(x)
    geno_grad = std::exp(log_like_cgcg_geno + LOG_PRIORS[0]) * reads_grad;

    // g'(x) (of the first derivative) -> z(x)*y'(x) + z'(x) * y(x)
    geno_grad_der2 = std::exp(log_like_cgcg_geno + LOG_PRIORS[0]) * reads_der2 + (geno_grad * reads_grad);
    log_like_genos = log_like_cgcg_geno + LOG_PRIORS[0];
    // summation of the likelihood of the remaining genotypes
    for (auto it=site->remaining_dinucl_genotypes.begin(); it!=site->remaining_dinucl_genotypes.end(); it++){
       log_like_genos = oplusnatl(log_like_genos, *it);
    }

    // f(x) first derivative is 1/below
    like_genos = std::exp(log_like_genos);
    // f'(x) of the first derivative
    like_genos_der2 = (- 1.0/(std::pow(like_genos,2))) * geno_grad;

    //             f'(x)           g(x)           f(x)              g'(x)
    total_d2 += like_genos_der2 * geno_grad + (1.0/like_genos) * geno_grad_der2;
  }
  return total_d2;
}


void setup_mma(nlopt::opt &opt){
  std::vector<double> lb(1, SMALLTOLERANCE), up(1, 1-SMALLTOLERANCE);

  opt.set_maxeval(1000);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(up);
  opt.set_xtol_abs(0);
  opt.set_ftol_abs(1e-15);
  opt.set_xtol_rel(0);
  opt.set_ftol_rel(0);

}

void print_header(std::ofstream &f, const std::string & l){
  f<<l;
}

boot_res_s calc_boots_f(general_settings & settings,  per_mle_run &mle_run,  F_void &void_stuff){
  nlopt::opt opt(nlopt::LD_MMA, 1);
  setup_mma(opt);

  std::vector<std::unique_ptr<Site_s>> boots;
  std::uniform_int_distribution<int> dis(0, mle_run.positions.size()-1);
  boots.reserve(mle_run.positions.size());
  std::vector<double> boot_res;

  for (size_t nb=0; nb<settings.nboots; nb++){
    for (size_t ns=0; ns<mle_run.positions.size(); ns++){
      size_t site_idx = dis(rn_generator);
      boots.emplace_back(std::make_unique<Site_s>(*void_stuff.to_include[site_idx]));
    }

    F_void void_stuff_boot(boots);
    std::vector<double> param_boot (1, SMALLTOLERANCE);
    double minf_boot;

    opt.set_max_objective(objective_func_F, &void_stuff_boot);
    nlopt::result result_boot = opt.optimize(param_boot, minf_boot);
    // std::cerr << nb << " " << minf_boot << " " << param_boot[0] << '\n';
    boot_res.push_back(param_boot[0]);
    boots.clear();
  }
  // calc summary statistics
  boot_res_s a;
  a.mean_boot = std::accumulate(boot_res.begin(), boot_res.end(), 0.0) / boot_res.size();
  for(auto &val: boot_res)
    a.var_boot += std::pow(val-a.mean_boot,2);
  a.var_boot /= boot_res.size();
  a.sd_boot = std::sqrt(a.var_boot);
  return (a);
}

void check_param_space(double &minimum_param, double &maximum_param){
  if(minimum_param<0)
    minimum_param=0;

  if(maximum_param>1)
    maximum_param = 1;
}


void run_mle(general_settings & settings,
             std::string & chrom,
             std::vector<Site_s> & cpg_data) {

  std::string filename = settings.outbase + "." + chrom + ".F";
  mtx.lock();
  settings.buffer += "\t-> Dumping MLE of F to " + filename + '\n';
  print_log(settings);
  mtx.unlock();
  std::ofstream f (filename.c_str());
  checkfilehandle<std::ofstream>(f, filename);

  if(settings.nboots)
    print_header(f, "center n_cpgs n_obs distance ll f sd_boot iter opt_return_code");
  else
    print_header(f, "center n_cpgs n_obs distance ll f f_95_conf iter opt_return_code");

  if(VERBOSE)
    print_header(f, " incl_pos\n");
  else
    print_header(f, "\n");


  std::set<size_t> last_positions_set;

  double last_minf;
  std::vector<double>  last_param;
  nlopt::result last_result;
  double last_error;
  double last_minimum_param, last_maximum_param;
  size_t last_iterations;

  // nlopt::opt opt(nlopt::LN_BOBYQA, 1);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  setup_mma(opt);

  for (size_t curr_idx = 0; curr_idx < cpg_data.size(); curr_idx++) {
    per_mle_run mle_run(cpg_data[curr_idx]);


    size_t pos_to_left = (curr_idx > 0) ? curr_idx - 1 : 0;
    size_t pos_to_right = curr_idx + 1;
    bool can_go_left = false, can_go_right = false;
    // keep adding sites to mle_run until we fulfill the requirements
    while (mle_run.max_pos-mle_run.min_pos < settings.windowsize && mle_run.n_cpgs < settings.max_cpgs) {
      can_go_right = false;
      can_go_left = false;

      if (pos_to_right < cpg_data.size() - 1) {
        can_go_right = true;
      }

      if (pos_to_left >= 1) {
        can_go_left = true;
      }

      if(std::max(mle_run.max_pos, mle_run.min_pos) - cpg_data[pos_to_left].pos > settings.windowsize){
        can_go_left = false;
      }

      if(cpg_data[pos_to_right].pos - std::min(mle_run.max_pos, mle_run.min_pos) > settings.windowsize){
        can_go_right = false;
      }

      if (!can_go_right && !can_go_left) {
        // break ;; cannot extend anymore
        break;

      } else if (can_go_right && !can_go_left) {
        // take the one to the right
        mle_run.update_right(cpg_data[pos_to_right]);
        //update_mle_run(mle_run, cpg_data[pos_to_right], pos_to_right);
        pos_to_right++;
      } else if (!can_go_right && can_go_left) {
        // take the one to the left
        mle_run.update_left(cpg_data[pos_to_left]);
        pos_to_left--;
      } else if (can_go_right && can_go_left) {
        // take the closest
        if (cpg_data[pos_to_right].pos - mle_run.curr_pos <=
            mle_run.curr_pos - cpg_data[pos_to_left].pos) {
          // take the one to the right
          mle_run.update_right(cpg_data[pos_to_right]);
          pos_to_right++;
        } else {
          // take the one to the left
          mle_run.update_left(cpg_data[pos_to_left]);
          pos_to_left--;
        }
      } else {
        // we are screwed. should not be possible
      }

    } // end while loop

    // here we can check if the results are identical to last time, then we dont need to do anything besides printing.
    // see below for the old code where we re calculate every time
    bool same = false;
    if(!last_positions_set.empty() && mle_run.positions.size()==last_positions_set.size()){
      same=true;
      for(auto &val: mle_run.positions){
        if(last_positions_set.count(val)!=1){
          same=false;
          break;
        }
      }
    }

    if(!same){
      last_positions_set.clear();
      last_positions_set.insert(mle_run.positions.begin(), mle_run.positions.end());
    }

    double minf;
    std::vector<double> param (1, SMALLTOLERANCE);
    nlopt::result result;
    double second_der, error;
    double minimum_param, maximum_param;
    size_t iterations;

    if(! same){
      F_void void_stuff(&settings, mle_run.to_include);
      opt.set_max_objective(objective_func_F, &void_stuff);
      result = opt.optimize(param, minf);

      if(settings.nboots){
        boot_res_s b_res = calc_boots_f(settings, mle_run, void_stuff);
        minimum_param = param[0]-b_res.sd_boot;
        maximum_param = param[0]+b_res.sd_boot;
        check_param_space(minimum_param, maximum_param);

      } else {
        second_der = objective_func_F_second_deriv(param[0], &void_stuff);
        error = 1.96/std::sqrt(-second_der);
        minimum_param = param[0]-error;
        maximum_param = param[0]+error;
        check_param_space(minimum_param, maximum_param);
      }

      iterations=void_stuff.iteration;

      // assign params
      last_minf = minf;
      last_param = param;
      last_result = result;
      last_minimum_param = minimum_param;
      last_maximum_param = maximum_param;
      last_iterations = iterations;
    } else {
      minf = last_minf;
      param = last_param;
      result = last_result;
      minimum_param = last_minimum_param;
      maximum_param = last_maximum_param;
      iterations = last_iterations;
    }

    // // keep this
    // double minf;
    // F_void void_stuff(settings, &cpg_data, &mle_run);
    // std::vector<double> param (1, 0.5);
    // nlopt::result result;
    // double second_der, error;
    // size_t iterations;
    // opt.set_max_objective(objective_func_F, &void_stuff);
    // result = opt.optimize(param, minf);
    // second_der = objective_func_F_second_deriv(param[0], cpg_data, mle_run);
    // error = 1.96/std::sqrt(-second_der);
    // iterations=void_stuff.iteration;
    // // to this


    // // temp
    // std::vector<double> dummy;
    // for (double i=0; i<=1; i+=0.01){
    //   std::vector<double> x(i);
    //   std::cout << i << " " << std::setprecision(10) << objective_func_F(param, dummy, &void_stuff) << std::endl;
    // }
    // // temp end

    // // probability of dinucleotide genotypes for the site in the center.
    // double first_dinucl_genotype = std::exp(loglikelihood_after_mle(param, &cpg_data[curr_idx]));
    // double sum_exp = first_dinucl_genotype;
    // for (auto & val : cpg_data[curr_idx].remaining_dinucl_genotypes){
    //   sum_exp += std::exp(val);
    // }
    // end

#if 0
    std::cerr << "done " << result << std::endl;
#endif

    f << chrom
      << ":" << mle_run.curr_pos
      << " " << mle_run.n_cpgs
      << " " << mle_run.total_depth
      << " " << mle_run.max_pos - mle_run.min_pos
      << " " << minf
      << " " << param[0]
      << " " << minimum_param  << "," <<  maximum_param
      << " " << iterations
      << " " << result;
    if(VERBOSE){
      f << " " << mle_run.positions[0];
      for (auto it=mle_run.positions.begin()+1; it!=mle_run.positions.end(); it++){
        f << ","<< *it;
      }
    }
      // // printing probability of dinucleotide genotypes for the site in the center.
      // f << " dinucl_genos: " << first_dinucl_genotype/sum_exp ;
      // for (auto & val : cpg_data[curr_idx].remaining_dinucl_genotypes){
      //        f << ","<< std::exp(val)/sum_exp ;
      // }
      // end

      f << '\n';

    // // // checking that the derivative is correct
    // double tmp_ll;
    // for (double fval=0.0; fval<=1; fval+=0.01){
    //   std::vector<double> f_vec_temp(1,fval), grad_vec_temp(1,0);
    //   tmp_ll = objective_func_F(f_vec_temp, grad_vec_temp, &void_stuff);
    //   f << "f: " << fval << " like: " << tmp_ll  << " " << std::exp(tmp_ll) << " derivative: "  << grad_vec_temp[0] << '\n';
    // }
    // f << std::flush;
    // f.close();

    // std::cerr << "BREAKING";
    // std::cerr << std::flush;
    // exit(EXIT_SUCCESS);


  }
  f << std::flush;
  f.close();
}


void run_mle_bed(general_settings & settings,
                 std::string & chrom,
                 std::vector<Site_s> & cpg_data) {
  std::vector<std::pair<size_t, size_t>> bed_coord = parse_bed_file(settings.bed_f,  chrom);
  std::string filename = settings.outbase + "." + chrom + ".BED.F";
  mtx.lock();
  settings.buffer += "\t-> Dumping MLE of F to " + filename + '\n';
  print_log(settings);
  mtx.unlock();
  std::ofstream f (filename.c_str());
  checkfilehandle<std::ofstream>(f, filename);

  // print header
  print_header(f, "bed n_cpgs n_obs ll f f_95_conf iter opt_return_code");
  if(VERBOSE)
    print_header(f, " incl_pos\n");
  else
    print_header(f, "\n");


  double last_minf;
  std::vector<double>  last_param;
  nlopt::result last_result;
  double last_error;
  size_t last_iterations;

  // nlopt::opt opt(nlopt::LN_BOBYQA, 1);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  setup_mma(opt);

  for (const auto & bed : bed_coord){
    std::vector<size_t> sites_to_include;
    for (size_t curr_idx = 0; curr_idx < cpg_data.size(); curr_idx++) {

      if (cpg_data[curr_idx].pos >= bed.first && cpg_data[curr_idx].pos < bed.second){

        sites_to_include.push_back(curr_idx);
      }
      if(cpg_data[curr_idx].pos >= bed.second){
        break;
      }

    }
    if( sites_to_include.size()==0){
      f << chrom << ":" << bed.first<<"-"<<bed.second << " " << "NO_SITES_IN_THE_REGION NOT_USED" << '\n';
      continue;
    }

    // it is stupid to run over data again, but too tired to make new structs for the BED setup.
    per_mle_run mle_run(cpg_data[sites_to_include[0]]);
    for (auto it=sites_to_include.begin()+1; it!=sites_to_include.end();it++){
      mle_run.stats_update(cpg_data[*it]);
    }

    if(mle_run.total_depth == 0){
      f << chrom << ":" << bed.first<<"-"<<bed.second << " " << "NO_DATA_IN_THE_REGION NOT_USED"  << '\n';
      continue;
    }

    double minf;
    F_void void_stuff(&settings, mle_run.to_include);
    std::vector<double> param (1, SMALLTOLERANCE);
    opt.set_max_objective(objective_func_F, &void_stuff);
    nlopt::result result = opt.optimize(param, minf);
    double second_der = objective_func_F_second_deriv(param[0], &void_stuff);
    double error = 1.96/std::sqrt(-second_der);


    f << chrom << ":" << bed.first << "-" << bed.second
      << " " << mle_run.n_cpgs
      << " " << mle_run.total_depth
      << " " << minf
      << " " << param[0]
      << " " << param[0]-error  << "," <<  param[0]+error
      << " " << void_stuff.iteration
      << " " << result;
    if(VERBOSE){
      f << " " << mle_run.positions[0];
      for (auto it=mle_run.positions.begin()+1; it!=mle_run.positions.end(); it++){
        f << ","<< *it;
      }
    }
    f << '\n';
  }
  f << std::flush;
  f.close();
}

void update_dinucl_priors(general_settings & settings){
  std::vector<double> res;
  std::stringstream ss( settings.priors_str );
  double val, sum=0;
  while( ss.good() ){
    std::string substr;
    getline( ss, substr, ',' );
    val = std::stod(substr);
    sum+=val;
    res.push_back( val );
  }
  if (res.size() != LOG_PRIORS.size()){
    std::cerr << "\t -> ERROR: provide a string of seven priors e.g. -h 1,0,0,0,0,0,0" << '\n';
    exit(EXIT_FAILURE);
  }

  // for (const auto & val : LOG_PRIORS){
  //   std::cerr << val << '\n';
  // }

  for (size_t i=0; i<res.size(); i++){
    LOG_PRIORS[i] = std::log(res[i]/sum);
  }
}

void load_rg_from_file(general_settings & settings, rgs_info &rgs){
    std::ifstream f (settings.readgroups_f.c_str());
    checkfilehandle<std::ifstream>(f, settings.readgroups_f);
    std::string row, rg, cycle_string;
    size_t cycle;
    std::stringstream ss;
    rgs.rg_split = true;
    while(getline(f, row)){
      ss.str(row);
      ss >> rg >> cycle_string;
      if(cycle_string.empty()){
        cycle=std::numeric_limits<size_t>::max();
        settings.buffer += "\t-> RG: " + rg + '\n';
      } else {
        cycle=std::stoi(cycle_string);
        settings.buffer += "\t-> RG: " + rg + ". Sequencing cycles: " + cycle_string + '\n';
      }
      rgs.rgs.push_back(rg);
      rgs.cycles.push_back(cycle);
      ss.clear();
      rg.clear();
      cycle_string.clear();
    }
    f.close();
    rgs.n = rgs.rgs.size();
}



struct read_group_info {
  std::string rgname;
  size_t rgname_idx;
};


// This might be faster:
// https://github.com/samtools/samtools/blob/72d140b590cbacc975e96bf40f3db6e6370a5cbe/bam.c#L77
int find_rg_idx(bam1_t *rd, rgs_info &rgs, read_group_info &read_rg){

  uint8_t * rgptr = bam_aux_get(rd, "RG");
  if (rgptr == NULL) {
    // in the case norg exists. read will be trashed
    // UNKWOWN
    std::string rname = std::string(bam_get_qname(rd));
    std::cerr << "[No_found_RG] :: " << rname << '\n';
    return -1;
  } else {
    // moves past the Z by +1. to get RG id
    read_rg.rgname = std::string((const char *)(rgptr + 1));
    auto match = std::find(rgs.rgs.begin(), rgs.rgs.end(), read_rg.rgname);
    if (match != rgs.rgs.end()) {
      read_rg.rgname_idx = std::distance(rgs.rgs.begin(), match);
      return 0;
    } else { // in the case it cannot find the rg. The read will be trashed
      return -1;
    }
  }
}

struct init_bam_s {
  samFile *in;
  bam_hdr_t *header;
  hts_idx_t *idx;
  hts_itr_t *iter;

  ~init_bam_s(){
    if(VERBOSE)
      std::cerr << "\t-> Freeing bams pointers" << '\n';

    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(in);
  }
};

std::unique_ptr<init_bam_s> my_init_bam(general_settings &settings, std::string &chrom){
  std::unique_ptr<init_bam_s> abc = std::make_unique<init_bam_s>();
  // open BAM for reading
  abc->in = sam_open(settings.bam_fn.c_str(), "r");
  if (abc->in == NULL) {
    std::cerr << "Unable to open BAM/SAM file: " << settings.bam_fn << '\n';
    exit(EXIT_FAILURE);
  }

  // Get the header
  abc->header = sam_hdr_read(abc->in);
  if (abc->header == NULL) {
    sam_close(abc->in);
    std::cerr << "Unable to open BAM header: " << settings.bam_fn << '\n';
    exit(EXIT_FAILURE);
  }

  // Load the index
  abc->idx = sam_index_load(abc->in, settings.bam_fn.c_str());
  if (abc->idx == NULL) {
    std::cerr
      << "Unable to open BAM/SAM index. Make sure alignments are indexed."
      << '\n';
    exit(EXIT_FAILURE);
  }

  // jump to chromosome
  abc->iter = sam_itr_querys(abc->idx, abc->header, chrom.c_str());
  if (abc->iter == NULL) {
    std::cerr
      << "Unable to jump to location."
      << '\n';
    exit(EXIT_FAILURE);
  }
  return abc;
}

struct init_ref_s {
  faidx_t *fai;
  char *ref;
  size_t seq_len;

  ~init_ref_s(){
    if(VERBOSE)
      std::cerr << "\t-> Freeing ref pointers stuff" << '\n';


    free(ref);
    free(fai);
  }
};

std::unique_ptr<init_ref_s> my_init_ref(general_settings &settings,
                       std::string &chrom){
  std::unique_ptr<init_ref_s> abc = std::make_unique<init_ref_s>();
  abc->fai = ref_init(settings.reference_fn);
  abc->ref = fetch_chrom(abc->fai, chrom);
  abc->seq_len = faidx_seq_len(abc->fai, chrom.c_str());
  return abc;
}

bool read_overlap_cpg(bam1_t *rd, const v_un_ch &cpg_bool){
  bool skipread = true;
  size_t startpos = rd->core.pos, endpos = bam_endpos(rd);

  for (size_t i=startpos; i<=endpos; i++){
    if(cpg_bool[i]==1){
      skipread=false;
      break;
    }
  }
  return skipread;
}

bool check_read_qual(const general_settings &settings, bam1_t *rd){
  bool remove=false;
  if (rd->core.l_qseq < settings.minreadlength ||
      rd->core.qual < settings.minmapQ ||
      (rd->core.flag & settings.flags_off) != 0 ||
      rd->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY))
    remove=true;

  return remove;
}

void mask_sites(general_settings &settings,
                std::string &chrom,
                std::unique_ptr<init_ref_s> &ref_s){

  if (!settings.exclude_bed_fn.empty()) {
    settings.buffer += "\t-> Masking genomic BED regions (-e) " + settings.exclude_bed_fn + '\n';
    filter_ref_bed(chrom, settings.exclude_bed_fn, ref_s->ref);
  }

  if (!settings.exclude_sites_fn.empty()) {
    settings.buffer += "\t-> Masking genomic sites (-E) " + settings.exclude_sites_fn + '\n';
    filter_ref_sites(chrom, settings.exclude_sites_fn, ref_s->ref);
  }
  print_log(settings);
}

struct read_stats {
  size_t counter=0, trashed=0, reads_skipped=0;
};


int check_and_align_read(general_settings &settings,
                         bam1_t *rd,
                         const v_un_ch &cpg_bool,
                         char *ref,
                         rgs_info &rgs,
                         read_stats &rs,
                         read_group_info &read_rg,
                         alignment_data &d){

    if (rs.counter && rs.counter % READS == 0) {
      std::cerr << "\t-> " << rs.counter << " reads processed and " << rs.trashed << " discarded. " << '\r';
    }
    rs.counter++;

    // check that the read overlaps a CpG
    if(read_overlap_cpg(rd, cpg_bool)){
      rs.reads_skipped++;
      return -1;
    }

    if(check_read_qual(settings, rd)){
      rs.trashed++;
      return -1;
    }

    if (rgs.rg_split) {
      int r = find_rg_idx(rd, rgs, read_rg);
      if(r<0){
        rs.trashed++;
        return -1;
      }
    } else {
      read_rg.rgname=ALL_RG;
      read_rg.rgname_idx=ALL_DAMMET_RG_IDX;
    }
    d=align_read(rd, ref);
    return 0;
}


void parse_reads_per_chrom_deamrates(general_settings & settings,
                                     std::string & chrom,
                                     std::vector<std::vector<int>> &tm,
                                     std::vector<uni_ptr_obs> &cpg_data,
                                     std::vector<uni_ptr_obs> &nocpg_data,
                                     rgs_info &rgs,
                                     std::vector<my_cov_rg> &cov_rg){

  std::unique_ptr<init_bam_s> bam_s = my_init_bam(settings, chrom);
  std::unique_ptr<init_ref_s> ref_s = my_init_ref(settings, chrom);
  mask_sites(settings, chrom, ref_s);
  const keeplist_map cpg_map = get_cpg_chrom_pos(ref_s->ref, ref_s->seq_len);
  int ncpgs=0;

  const v_un_ch cpg_bool = get_cpg_chrom(ref_s->ref, ref_s->seq_len, ncpgs);


  read_stats rs;
  read_group_info read_rg;
  alignment_data d;
  bam1_t *rd = bam_init1();

  time_t start_time_load_data,end_time_load_data;
  time(&start_time_load_data);

  while (sam_itr_next(bam_s->in, bam_s->iter, rd) >= 0) {
    int r = check_and_align_read(settings, rd, cpg_bool, ref_s->ref, rgs, rs, read_rg, d);
    if(r<0)
      continue;
    add_aligned_data(settings, d,
                     cpg_data[read_rg.rgname_idx],
                     nocpg_data[read_rg.rgname_idx],
                     tm[read_rg.rgname_idx],
                     cov_rg[read_rg.rgname_idx],
                     rgs.cycles[read_rg.rgname_idx]);
  }
  time(&end_time_load_data);
  mtx.lock();
  settings.buffer += "\t-> " + std::to_string(ncpgs) + " CpG's in chrom: " + chrom + '\n';
  settings.buffer += "\t-> Chrom: " + chrom + ". Processed: " +
    std::to_string(rs.counter)  + ". Reads filtered: " + std::to_string(rs.trashed) +
    ". Reads skipped (nocpg overlap): " + std::to_string(rs.reads_skipped) +
    ". Loaded in " + std::to_string(difftime(end_time_load_data, start_time_load_data)) + " seconds." + '\n';
  print_log(settings);
  mtx.unlock();
  // print_log(settings);
  bam_destroy1(rd);
}

void parse_deamrates_wrapper(general_settings & settings, job_deamrates &jb){
  for(auto & c : jb.chroms){
    parse_reads_per_chrom_deamrates(settings, c, jb.tm,
                                    jb.cpg_data, jb.nocpg_data,
                                    jb.rgs, jb.cov_rg);
  }
}


std::vector<Site_s> remove_cpg_wo_data(std::vector<Site_s> &data){
  std::vector<Site_s> res;
  for (auto &val: data){
    if(val.depth)
      res.push_back(val);
  }
  return(res);
}

void parse_reads_per_chrom_estF(general_settings & settings,
                                std::string & chrom,
                                std::vector<std::vector<double>> &param_deam_rgs,
                                rgs_info &rgs){


  std::unique_ptr<init_bam_s> bam_s = my_init_bam(settings, chrom);
  std::unique_ptr<init_ref_s> ref_s = my_init_ref(settings, chrom);
  mask_sites(settings, chrom, ref_s);

  const keeplist_map cpg_map = get_cpg_chrom_pos(ref_s->ref, ref_s->seq_len);
  int ncpgs=0;
  // const std::vector<int> cpg_bool = get_cpg_chrom_bool(ref, seq_len);
  const v_un_ch cpg_bool = get_cpg_chrom(ref_s->ref, ref_s->seq_len, ncpgs);

  read_stats rs;
  read_group_info read_rg;
  std::vector<my_cov_rg> cov_rg;
  cov_rg.resize(rgs.n);

  alignment_data d;
  std::vector<Site_s> cpg_data;
  cpg_data.resize(ncpgs);
  for (auto &m: cpg_map){
    cpg_data[m.second].pos = m.first;
  }

  bam1_t *rd = bam_init1();

  time_t start_time_load_data,end_time_load_data;
  time(&start_time_load_data);

  while (sam_itr_next(bam_s->in, bam_s->iter, rd) >= 0) {
    int r = check_and_align_read(settings, rd, cpg_bool, ref_s->ref, rgs, rs, read_rg, d);
    if(r<0)
      continue;

    // do the work here
    calc_M_noM(settings,
               d,
               cpg_map,
               rgs.cycles[read_rg.rgname_idx],
               cov_rg[read_rg.rgname_idx],
               param_deam_rgs[read_rg.rgname_idx],
               cpg_data);

  }
  time(&end_time_load_data);

  mtx.lock();
  settings.buffer += "\t-> " + std::to_string(ncpgs) + " CpG's in chrom: " + chrom + '\n';
  for (size_t i=0; i<rgs.n; i++){
     settings.buffer += "\t-> Total_Observations RG: " + rgs.rgs[i] + " CpG obs: " +
       std::to_string(cov_rg[i].cpg) + " CpGCoverage: " +
       std::to_string((double)cov_rg[i].cpg/(double)ncpgs) + '\n';
  }
  print_log(settings);
  mtx.unlock();
  // adding priors to all alt genotypes.
  for(auto & s: cpg_data){
    for (size_t dinucl_idx=1; dinucl_idx<SEVEN_DINUCL_GENOTYPES.size(); dinucl_idx++){
      s.remaining_dinucl_genotypes[dinucl_idx-1] += s.remaining_dinucl_genotypes[dinucl_idx-1] + LOG_PRIORS[dinucl_idx];
    }
  }


  if(settings.skip_empty_cpg){
    mtx.lock();
    cpg_data = remove_cpg_wo_data(cpg_data);
    settings.buffer = "\t-> Skip CpGs without data. CpGs remaining: " + std::to_string(cpg_data.size()) + '\n';
    print_log(settings);
    mtx.unlock();
  }


  if(settings.bed_f.empty()){
    run_mle(settings, chrom, cpg_data);
  } else {
    run_mle_bed(settings, chrom, cpg_data);
  }


  bam_destroy1(rd);
}

void parse_fest_wrapper(general_settings & settings, job_fest &jb){
  for(auto & c : jb.chroms){
    parse_reads_per_chrom_estF(settings, c, jb.param_deam_rgs, jb.rgs);
  }
}



void compare_chroms(const general_settings &settings,
                    std::vector<std::string> &chroms,
                    const std::string &filename){
    for(auto &c : settings.chrom){
      auto match = std::find(chroms.begin(), chroms.end(), c);
      if (match == chroms.end()) {
        std::cerr << "\nChrom '" << c << "' not found in '" << filename << '\n';
        std::cerr << "EXITING" << '\n';
        exit(EXIT_FAILURE);
      }
    }
}

void chrom_in_fai(const general_settings &settings){
  std::string filename = settings.reference_fn + ".fai";
  std::ifstream f (filename.c_str());
  checkfilehandle<std::ifstream>(f, filename);
  std::string row, chrom;
  std::vector<std::string> fai_chroms;
  while(getline(f,  row)){
    std::stringstream ss(row);
    ss >> chrom;
    fai_chroms.push_back(chrom);
  }
  compare_chroms(settings, fai_chroms, filename);

  f.close();
}

void chrom_in_bam(const general_settings &settings){
  // open BAM for reading
  samFile *in = sam_open(settings.bam_fn.c_str(), "r");
  if (in == NULL) {
    std::cerr << "Unable to open BAM/SAM file: " << bam << '\n';
    exit(EXIT_FAILURE);
  }

  // Get the header
  bam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL) {
    sam_close(in);
    std::cerr << "Unable to open BAM header: " << bam << '\n';
    exit(EXIT_FAILURE);
  }

  std::vector<std::string> ref_names;
  for (int i=0; i<header->n_targets; i++){
    ref_names.push_back(std::string(header->target_name[i]));
  }

  bam_hdr_destroy(header);
  sam_close(in);
  compare_chroms(settings, ref_names, settings.bam_fn);
}

void merge_threads_deamrates(std::vector<job_deamrates> &jobs, job_deamrates & res){
  for(auto &j:jobs){
    for (size_t i=0; i<res.rgs.n; i++){

      // tally counts
      for (size_t ii=0; ii<j.tm[i].size(); ii++)
        res.tm[i][ii] += j.tm[i][ii];

      // coverage
      res.cov_rg[i].nocpg +=  j.cov_rg[i].nocpg;
      res.cov_rg[i].cpg +=  j.cov_rg[i].cpg;

      res.cpg_data[i].insert(res.cpg_data[i].end(),
                  std::make_move_iterator(j.cpg_data[i].begin()),
                  std::make_move_iterator(j.cpg_data[i].end()));
      res.nocpg_data[i].insert(res.nocpg_data[i].end(),
                  std::make_move_iterator(j.nocpg_data[i].begin()),
                  std::make_move_iterator(j.nocpg_data[i].end()));
    }
  }
}

job_deamrates compute_multithreading_deamrates(general_settings & settings, rgs_info &rgs){
  size_t block = settings.nthreads==1?settings.chrom.size(): settings.chrom.size() / settings.nthreads;
  std::vector<job_deamrates> jobs;
  size_t curr_chrom_idx=0;
  for(size_t i=0; i<settings.nthreads; i++){
    jobs.push_back(job_deamrates(rgs));
    for (size_t iii=0; iii<rgs.n; iii++)
      jobs[i].tm[iii] = init_tallymat<int>(settings.max_pos_to_end);

    int chrom_idx =i*block;
    for(size_t ii=chrom_idx; ii<(chrom_idx+block) && ii<settings.chrom.size(); ii++){
      jobs[i].chroms.push_back(settings.chrom[ii]);
      curr_chrom_idx++;
    }
  }
  size_t n_extra = settings.chrom.size() % settings.nthreads;
  if(n_extra){
    std::cerr << "\t-> Adding remaining chromosomes " << n_extra << '\n';

    while (n_extra){

      jobs[jobs.size()-n_extra].chroms.push_back(settings.chrom[curr_chrom_idx]);
      curr_chrom_idx++;
      n_extra--;
    }
  }
#if 0
  std::cerr << settings.chrom.size() << std::endl;
  std::cerr << "entering" << '\n';
  for(size_t i=0; i<settings.nthreads; i++){
    int chrom_idx =i*block;
    for(size_t ii=0; ii<jobs[i].chroms.size(); ii++){
      std::cerr << i << " " << ii <<  " " << jobs[i].chroms[ii] << '\n';;
    }
  }
  std::cerr << std::flush;
  exit(0);
#endif



  std::vector<std::thread> threads;
  for(size_t i=0; i<settings.nthreads; i++){
    // std::thread a threadobj(parse_deamrates_wrapper, jobs[i]);
    // threads.push_back(std::move(a));
    threads.push_back(std::thread(parse_deamrates_wrapper, std::ref(settings), std::ref(jobs[i])));
  }

  for(size_t i=0; i<settings.nthreads; i++){
    if(threads[i].joinable())
      threads[i].join();

  }

  job_deamrates res(rgs);
  for (size_t iii=0; iii<rgs.n; iii++)
      res.tm[iii] = init_tallymat<int>(settings.max_pos_to_end);

  merge_threads_deamrates(jobs, res);
  return(res);

}


void estdeam(general_settings & settings, rgs_info &rgs) {
  print_log(settings);
  time_t mu_st, mu_et;
  time(&mu_st);
  job_deamrates merged = compute_multithreading_deamrates(settings, rgs);
  time(&mu_et);
  double mu_tg = difftime(mu_et, mu_st);

  std::cerr << "\t-> Done loading " << settings.chrom.size() << " chromosomes in " << mu_tg << " seconds." << '\n';
  for (size_t i=0; i<rgs.n; i++){
    settings.buffer += "\t-> Total_Observations RG: "+ rgs.rgs[i] +
      " CpG/cnonCpG: " + std::to_string(merged.cov_rg[i].cpg) + " " + std::to_string(merged.cov_rg[i].nocpg) + '\n';
    print_log(settings);
  }

  // exit(0);
  std::vector<std::vector<double>> tmf;
  tmf.resize(rgs.n);
  // dumping counts to file
  for (size_t i=0; i<rgs.n; i++){

    settings.buffer += "\t-> Dumping count file: " + settings.outbase + "." + merged.rgs.rgs[i] +  ".tallycounts" + '\n';
    print_log(settings);
    dump_count_file(settings, merged.tm[i], rgs.rgs[i]);
    tmf[i] = read_count_file(settings, rgs.rgs[i]);
  }

  std::vector<std::vector<double>> param_deam_rgs;
  param_deam_rgs.resize(rgs.n);
  settings.args_stream << std::flush;

  for (size_t i=0; i<rgs.n; i++){
    settings.buffer += "\t-> Starting Optim of deamination rates. RG: " + rgs.rgs[i] + '\n';
    print_log(settings);

#if 0
    print_data(merged.cpg_data[i], merged.nocpg_data[i]);
#endif
    run_deamrates_optim(settings, tmf[i], merged.cpg_data[i], merged.nocpg_data[i], rgs.rgs[i]);
    settings.buffer += "\t-> Dumping deamination parameters to " + settings.outbase+"."+rgs.rgs[i]+".deamrates" + '\n';
    print_log(settings);
  }
  settings.args_stream << std::flush;
}

void estF(general_settings & settings, rgs_info &rgs){

  // load deamination rates from deamrates
  std::vector<std::vector<double>> param_deam_rgs;
  param_deam_rgs.resize(rgs.n);
  for (size_t i=0; i<rgs.n; i++){
    if((!settings.deamrates_filename.empty()) && check_file_exists(settings.deamrates_filename)){
      settings.buffer += "\t-> Loading deamination rates from " + settings.deamrates_filename + '\n';
      print_log(settings);
      std::cerr << "\t-> Make sure that the file contains the same number of pos to include. DamMet does not check that" << '\n';
      param_deam_rgs[i] = load_deamrates_f(settings);
    } else if (check_file_exists(settings.outbase+"."+rgs.rgs[i]+".deamrates")){
      settings.buffer += "\t-> Loading deamination rates from " + settings.outbase+"."+rgs.rgs[i]+".deamrates" + '\n';
      print_log(settings);
      param_deam_rgs[i] = load_deamrates(settings, rgs.rgs[i]);
    } else{
      std::cerr << "\nCannot find deamination rates. EXITING" << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  print_log(settings);

  // update priors if provided
  if(!settings.priors_str.empty()){
    update_dinucl_priors(settings);
  }

  settings.buffer = "\t-> Dinucl genotype priors (-h): ";

  for (const auto & val : LOG_PRIORS){
    settings.buffer += std::to_string(val)+',';
  }
  settings.buffer += '\n';
  print_log(settings);



  size_t block = settings.nthreads==1?settings.chrom.size(): settings.chrom.size() / settings.nthreads;
  std::vector<job_fest> jobs;
  size_t curr_chrom_idx=0;
  for(size_t i=0; i<settings.nthreads; i++){
    jobs.push_back(job_fest(rgs, param_deam_rgs));

    int chrom_idx =i*block;
    for(size_t ii=chrom_idx; ii<(chrom_idx+block) && ii<settings.chrom.size(); ii++){
      jobs[i].chroms.push_back(settings.chrom[ii]);
      curr_chrom_idx++;
    }
  }
  size_t n_extra = settings.chrom.size() % settings.nthreads;
  if(n_extra){
    std::cerr << "\t-> Adding remaining chromosomes " << n_extra << '\n';

    while (n_extra){
      jobs[jobs.size()-n_extra].chroms.push_back(settings.chrom[curr_chrom_idx]);
      curr_chrom_idx++;
      n_extra--;
    }
  }

  std::vector<std::thread> threads;
  for(size_t i=0; i<settings.nthreads; i++){
    // std::thread a threadobj(parse_deamrates_wrapper, jobs[i]);
    // threads.push_back(std::move(a));
    threads.push_back(std::thread(parse_fest_wrapper, std::ref(settings), std::ref(jobs[i])));
  }

  for(auto &thread: threads){
    if(thread.joinable())
      thread.join();
  }
}


int main(int argc, char *argv[]) {
  time_t start_time, end_time;
  time(&start_time);
  general_settings settings;
  args_parser(argc, argv, settings);

  VERBOSE = settings.verbose;

  if(settings.seed != std::numeric_limits<size_t>::max()){
    rn_generator.seed(settings.seed);
  }

  if(settings.nthreads > settings.chrom.size()){
    std::cerr << "\t-> Less contigs than threads. Reducing nthreads." << '\n';
    settings.nthreads = settings.chrom.size();
  }

  if(settings.analysis=="estF"){
    if(check_estF_args(settings)){
      print_help();
      std::cerr << "Must specify either -N (N CpGs per window) AND/OR -W (max windowsize) OR -B bedfile" << '\n';
      std::cerr << "EXITING...." << '\n';
      exit(EXIT_FAILURE);
    }

    std::string stream_filename (settings.outbase+".estF.args");
    settings.args_stream.open(stream_filename.c_str());
    checkfilehandle<std::ofstream>(settings.args_stream, stream_filename);
  } else {
    std::string stream_filename (settings.outbase+".estDEAM.args");
    settings.args_stream.open(stream_filename.c_str());
    checkfilehandle<std::ofstream>(settings.args_stream, stream_filename);
  }

  settings.buffer += "\t-> nthreads: " + std::to_string(settings.nthreads) + '\n';
  print_log(settings);
  // print general info to args file.

  rgs_info rgs;

  if((!settings.readgroups_f.empty()) && check_file_exists(settings.readgroups_f)){
    load_rg_from_file(settings, rgs);
  } else {
    rgs.rg_split = false;
    rgs.rgs.push_back(ALL_RG);
    settings.buffer += "\t-> Merging all reads into a single read group named: " + ALL_RG;
    if(settings.cycles!=std::numeric_limits<size_t>::max()){
      rgs.cycles.push_back(settings.cycles);
      settings.buffer += ". Sequencing cycles: " + std::to_string(settings.cycles);
    } else {
      rgs.cycles.push_back(std::numeric_limits<size_t>::max());
    }
    settings.buffer += '\n';
    rgs.n = 1;
  }

  std::cerr << "\t-> DamMet started: " << ctime(&start_time);

  // check that chromosome is in .fai and .bam
  chrom_in_bam(settings);
  chrom_in_fai(settings);

  // print options to args file
  print_log(settings, false);

  if(settings.analysis=="estF"){
    settings.buffer += "\t-> Estimating methylation levels (f)\n";
    print_log(settings);
    estF(settings, rgs);
  } else {
    size_t files_avail = 0;
    for (size_t i=0; i<rgs.n; i++){

      if((!settings.deamrates_filename.empty()) && check_file_exists(settings.deamrates_filename)){
        settings.buffer += "\t-> As deamination profiles are provided to -D run 'estF' to estimate methylation levels\n" ;
        print_log(settings);
        std::cerr << "\t-> Make sure that the file contains the same number of pos to include. DamMet does not check that" << '\n';
        files_avail++;
      } else if (check_file_exists(settings.outbase+"."+rgs.rgs[i]+".deamrates")){
        settings.buffer += "\t-> Deamination rates from " + settings.outbase+"."+rgs.rgs[i]+".deamrates are available" + '\n';
        print_log(settings);
        files_avail++;
      }
    }

    if(files_avail != rgs.n){
      settings.buffer += "\t-> Estimating Deamination rates (D)\n";
      print_log(settings);
      estdeam(settings, rgs);
    } else {
      std::cerr << "\t-> Deamination rates are already calculated. Please, delete the files if you want to recalculate the deamination rates.\n";
    }
  }


  time(&end_time);
  double time_gone = difftime(end_time, start_time);
  size_t minutes = time_gone / 60.0;
  size_t seconds = (int)time_gone % 60;
  settings.buffer += "\t-> Done in " + std::to_string(minutes)+ ":" + std::to_string(seconds) + " M:S." + '\n';
  std::cerr << "\t-> DamMet ended: " << ctime(&end_time);
  print_log(settings);
}
