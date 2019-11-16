#include "version.hpp"

#include "dammet_dam_params.hpp"
#include "nlopt.hpp"
// #include "/home/krishang/Desktop/test2/DamMet/nlopt-2.5.0/install/include/nlopt.h"
#include <set>
#include <ctime>
#include <iomanip> // setprecision


template <typename T>
void checkfilehandle(T &fh, std::string filename){
  if (! fh.is_open()){
    std::cerr << "Couldnt open file: " << filename << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
}

template <class T>
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

void print_log(general_settings & settings){
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

std::vector<double> PHRED_TO_PROB_CONVERTER = phred_to_double_converter();


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

// const std::vector<int> get_cpg_chrom_pos (const char * ref, const size_t & seq_len){
//   size_t counter = 0;
//   for (size_t i=0; i<seq_len-1; i++){
//     if (refToInt[(int)ref[i]] == 1 && refToInt[(int)ref[i+1]] == 2){
//       res.push_back(i);
//     }
//   }
//   return res;
// }

// const v_un_ch get_c_and_cpg_chrom(const char * ref, const size_t & seq_len){
//   v_un_ch res(seq_len, 0);
//   for (size_t i=0; i<seq_len-1; i++){
//     if (refToInt[(int)ref[i]] == 1){

//       if (refToInt[(int)ref[i+1]] == 0){
//         res[i] = 1;
//         res[i+1] = 1;
//       } else if (refToInt[(int)ref[i+1]] == 1){
//         res[i] = 2;
//         res[i+1] = 2;
//       } else if (refToInt[(int)ref[i+1]] == 2){
//         res[i] = 3;
//         res[i+1] = 3
//       } else if (refToInt[(int)ref[i+1]] == 3){

//       }
//     }
//   }
//   return res;
// }

const v_un_ch get_cpg_chrom(const char * ref, const size_t & seq_len){
  v_un_ch res(seq_len, 0);
  for (size_t i=0; i<seq_len-1; i++){
    if ((refToInt[(int)ref[i]] == 1)&&(refToInt[(int)ref[i+1]] == 2)){
      res[i]=1;
      res[i+1]=1;
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

// std::vector<int> idx_to_params(const int & idx){
//   std::div_t p = std::div(idx, READPOS_MULT);
//   std::div_t pr = std::div(p.rem, PRIME_MULT);
//   std::div_t strand = std::div(pr.rem, STRAND_MULT);
//   std::div_t B1 = std::div(strand.rem, B1_MULT);
//   std::div_t B2 = std::div(B1.rem, B2_MULT);
//   std::div_t B3 = std::div(B2.rem, B3_MULT);
//   return std::vector<int> {p.quot, pr.quot, strand.quot, B1.quot, B2.quot, B3.quot, B3.rem};
// }

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

bool nucleotide_missing(const int & r_base1, const int & r_base2, const int & s_base1, const int & s_base2){
  return (r_base1>=4 || r_base2>=4 || s_base1>=4 || s_base2>=4);
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

void add_aligned_data(general_settings &settings,
                      const alignment_data &d,
                      uni_ptr_obs &cpg_data,
                      uni_ptr_obs &nocpg_data,
                      std::vector<int> &tm,
                      my_cov_rg &cov_rg,
                      size_t &rg_idx,
                      size_t &ncycles,
                      size_t &exclude_CnonCpG) {

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
         cov_rg.nocpg[rg_idx] < cov_rg.cpg[rg_idx]){
        nocpg_data.emplace_back(std::make_unique<Obs>(prime,
                                                      d.strand,
                                                      pos,
                                                      rl,
                                                      bc,
                                                      d.t_qs[boi],
                                                      d.mapQ));
        cov_rg.nocpg[rg_idx]++;
      }

      if(cpg(d.t_ref[boi-1], d.t_ref[boi])){
        cpg_data.emplace_back(std::make_unique<Obs>(prime,
                                                    d.strand,
                                                    pos,
                                                    rl,
                                                    bc,
                                                    d.t_qs[boi],
                                                    d.mapQ));
        cov_rg.cpg[rg_idx]++;
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

      if(no_ref_cpg(d.t_ref[boi], d.t_ref[boi+1]) && cov_rg.nocpg[rg_idx] < cov_rg.cpg[rg_idx]){
          // adding data to nocpg
          nocpg_data.emplace_back(std::make_unique<Obs>(prime, d.strand, pos, rl, bc, d.t_qs[boi], d.mapQ));
          cov_rg.nocpg[rg_idx]++;
      }

      if(cpg(d.t_ref[boi], d.t_ref[boi+1])){
        // adding data to cpg
        cpg_data.emplace_back(std::make_unique<Obs>(prime, d.strand, pos, rl, bc, d.t_qs[boi], d.mapQ));
        cov_rg.cpg[rg_idx]++;
      }

      tm[get_idx_tm(pos, prime, d.strand, d.t_ref[boi], d.t_ref[boi+1], d.t_seq[boi], d.t_seq[boi+1])]++;
    }
  }
}

size_t get_param_idx(const size_t & max_pos_to_end, const size_t & meth, const size_t & pos_to_end, const size_t & prime){
  return (meth * PRIMES * ( max_pos_to_end+2 ) + pos_to_end * PRIMES  + prime);
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

void update_mle_data(per_mle_run & mle_data, const pre_calc_per_site & data, size_t & idx){
  mle_data.idx_to_include.push_back(idx);
  mle_data.total_depth += data.depth;
  mle_data.positions.push_back(data.position);
  mle_data.n_cpgs++;
  if (mle_data.min_pos > data.position) {
    mle_data.min_pos = data.position;
  } else if (mle_data.max_pos < data.position) {
    mle_data.max_pos = data.position;
  } else {
    std::cerr << mle_data.curr_pos << " " << mle_data.min_pos << " " << mle_data.max_pos << " NOT GOOD POSITION " << data.position << '\n';
  }
}

double objective_func_F_haploid(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
  F_void *d = static_cast<F_void*>(my_func_data);
  d->iteration++;
  double res, ll=0;
  double noM, M;
  double mape;
  if(!grad.empty()){
    grad[0] = 0;
  }
  pre_calc_per_site * site;
  for (const auto idx : d->mle_data->idx_to_include){
    site = &d->data->at(idx);
    for (size_t i=0; i < site->depth; i++){
      noM =  site->pre_noM[i];
      M = site->pre_M[i];
      mape = site->maperrors[i];
      res = (1-mape) * ((1-x[0]) * noM + x[0] * M) + mape*BASE_FREQ_FLAT_PRIOR;
      ll += std::log(res);
      // wolfram: derivative  ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // slope += ((1-mape) * (M - noM)) / ((mape-1) * ((f-1) * noM - f * M) + mape*BASE_FREQ_FLAT_PRIOR);
      // beloved wolfram got it right, just modified a few things for readability
      if(!grad.empty()){
        grad[0] += ((1-mape) * (M - noM)) / res;
      }
      // https://www.symbolab.com/
      // wolfram: second derivative ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // second derivative:
      // second_der += (std::pow((1-mape),2) * std::pow((M - noM),2)) / std::pow((-M * (mape-1) * f + noM * (mape-1) * (1-f) + mape*BASE_FREQ_FLAT_PRIOR), 2);
    }
  }
  return ll;
}

double objective_func_F_second_deriv_haploid(const double & f, std::vector<pre_calc_per_site> & data, per_mle_run &mle_data){
  double slope=0;
  double nominator, denominator;
  double noM, M;
  double mape;
  pre_calc_per_site * site;
  for (const auto idx : mle_data.idx_to_include){
    site = &data.at(idx);
    for (size_t i=0; i < site->depth; i++){
      noM =  site->pre_noM[i];
      M = site->pre_M[i];
      mape = site->maperrors[i];
      nominator = std::pow((-mape + 1 ), 2) * std::pow((M - noM ), 2);
      denominator = std::pow((1-mape) * (noM * (1-f) + f * M) + mape*BASE_FREQ_FLAT_PRIOR, 2);
      slope += -nominator/denominator;
    }
  }
  return slope;
}


double objective_func_F(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
  F_void *d = static_cast<F_void*>(my_func_data);
  d->iteration++;
  double log_like_cgcg_geno, read_der;
  double reads_grad, geno_grad;
  double read_like, log_like_genos ;
  double noM, M, mape ;
  double log_like_window=0;
  if(!grad.empty()){
    grad[0] = 0;
  }
  pre_calc_per_site * site;
  for (const auto & idx : d->mle_data->idx_to_include){
    site = &d->data->at(idx);
    log_like_cgcg_geno = 0;
    reads_grad = 0;
    geno_grad = 0;
    log_like_genos = 0;
    for (size_t i=0; i < site->depth; i++){
      noM =  site->pre_noM[i];
      M = site->pre_M[i];
      mape = site->maperrors[i];
      read_der = (1-mape) * (M - noM);
      read_like = (1-mape) * ((1-x[0]) * noM + x[0] * M) + mape*DINUCL_FLAT_PRIOR;
      log_like_cgcg_geno += std::log(read_like);
      if(!grad.empty()){
        // d/df per read
        reads_grad += read_der/read_like;
      }
    // wolfram: derivative  ln((1-w)*((1-x) * k + (x * h)) + w*p)
    // beloved wolfram got it right, just modified a few things for readability
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

double loglikelihood_after_mle(const std::vector<double> &x, pre_calc_per_site * site){
  double log_like_cgcg_geno=0;
  double read_like ;
  double noM, M, mape ;
  for (size_t i=0; i < site->depth; i++){
    noM =  site->pre_noM[i];
    M = site->pre_M[i];
    mape = site->maperrors[i];
    read_like = (1-mape) * ((1-x[0]) * noM + x[0] * M) + mape*DINUCL_FLAT_PRIOR;
    log_like_cgcg_geno += std::log(read_like);
  }
  return log_like_cgcg_geno + LOG_PRIORS[0];
}

double objective_func_F_second_deriv(const double & f, std::vector<pre_calc_per_site> & data, per_mle_run &mle_data){
  double total_d2=0;
  // double nominator, denominator;
  double log_like_cgcg_geno;
  double reads_grad, geno_grad;
  double read_der, read_like, log_like_genos ;
  double noM, M;
  double mape;
  double numerator_read_der2, denominator_read_der2, reads_der2;
  double geno_grad_der2, like_genos, like_genos_der2;
  pre_calc_per_site * site;
  for (const auto & idx : mle_data.idx_to_include){
    site = &data.at(idx);
    log_like_cgcg_geno = 0;
    reads_grad=0;
    reads_der2=0;
    for (size_t i=0; i < site->depth; i++){
      noM =  site->pre_noM[i];
      M = site->pre_M[i];
      mape = site->maperrors[i];

      read_der = (1-mape) * (M - noM);
      read_like = (1-mape) * ((1-f) * noM + f * M) + mape*DINUCL_FLAT_PRIOR;
      reads_grad += read_der/read_like;

      numerator_read_der2 = std::pow((1-mape),2) * std::pow((M - noM),2);
      denominator_read_der2 = std::pow((1-mape) * (noM * (1-f) + f * M) + mape*DINUCL_FLAT_PRIOR, 2);
      reads_der2 += - (numerator_read_der2 / denominator_read_der2);
      // \frac{\partial }{\partial \:x^2}\left(ln\left(\left(1-w\right)\cdot \left(\left(1-x\right)\:\cdot \:n\:+\:\left(x\:\cdot \:M\right)\right)\:+\:w\cdot \:p\right)\:\:\right)

      log_like_cgcg_geno += std::log(read_like);
      // https://www.symbolab.com/
      // wolfram: second derivative ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // second derivative:
      // second_der += (std::pow((1-mape),2) * std::pow((M - noM),2)) / std::pow((-M * (mape-1) * f + noM * (mape-1) * (1-f) + mape*DINUCL_FLAT_PRIOR), 2);
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


void run_mle(general_settings & settings,
             std::vector<pre_calc_per_site> & pre_calc_data) {
  std::string filename = settings.outbase + ".F";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);

  std::set<size_t> last_positions_set;

  double last_minf;
  std::vector<double>  last_param;
  nlopt::result last_result;
  double last_error;
  size_t last_iterations;

  // nlopt::opt opt(nlopt::LN_BOBYQA, 1);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  std::vector<double> lb(1, SMALLTOLERANCE), up(1, 1-SMALLTOLERANCE);

  opt.set_maxeval(1000);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(up);
  opt.set_xtol_abs(0);
  opt.set_ftol_abs(1e-15);
  opt.set_xtol_rel(0);
  opt.set_ftol_rel(0);
  for (size_t curr_idx = 0; curr_idx < pre_calc_data.size(); curr_idx++) {
    per_mle_run mle_data(pre_calc_data[curr_idx], curr_idx);
    size_t pos_to_left = (curr_idx > 0) ? curr_idx - 1 : 0;
    size_t pos_to_right = curr_idx + 1;
    bool can_go_left = false, can_go_right = false;
    // keep adding sites to mle_data until we fulfill the requirements
    while (mle_data.max_pos-mle_data.min_pos < settings.windowsize && mle_data.n_cpgs < settings.max_cpgs) {
      can_go_right = false;
      can_go_left = false;

      if (pos_to_right < pre_calc_data.size() - 1) {
        can_go_right = true;
      }

      if (pos_to_left >= 1) {
        can_go_left = true;
      }

      if(std::max(mle_data.max_pos, mle_data.min_pos) - pre_calc_data[pos_to_left].position > settings.windowsize){
        can_go_left = false;
      }

      if(pre_calc_data[pos_to_right].position - std::min(mle_data.max_pos, mle_data.min_pos) > settings.windowsize){
        can_go_right = false;
      }

      if (!can_go_right && !can_go_left) {
        // break ;; cannot extend anymore
        break;

      } else if (can_go_right && !can_go_left) {
        // take the one to the right
        update_mle_data(mle_data, pre_calc_data[pos_to_right], pos_to_right);
        pos_to_right++;
      } else if (!can_go_right && can_go_left) {
        // take the one to the left
        update_mle_data(mle_data, pre_calc_data[pos_to_left], pos_to_left);
        pos_to_left--;
      } else if (can_go_right && can_go_left) {
        // take the closest
        if (pre_calc_data[pos_to_right].position - mle_data.curr_pos <=
            mle_data.curr_pos - pre_calc_data[pos_to_left].position) {
          // take the one to the right
          update_mle_data(mle_data, pre_calc_data[pos_to_right], pos_to_right);
          pos_to_right++;
        } else {
          // take the one to the left
          update_mle_data(mle_data, pre_calc_data[pos_to_left], pos_to_left);
          pos_to_left--;
        }
      } else {
        // we are screwed. should not be possible
      }

    } // end while loop

    // here we can check if the results are identical to last time, then we dont need to do anything besides printing.
    // see below for the old code where we re calculate every time
    bool same = false;
    if(!last_positions_set.empty() && mle_data.positions.size()==last_positions_set.size()){
      same=true;
      for(const size_t &  val: mle_data.positions){
        if(last_positions_set.count(val)!=1){
          same=false;
          break;
        }
      }
    }

    if(!same){
      last_positions_set.clear();
      last_positions_set.insert(mle_data.positions.begin(), mle_data.positions.end());
    }

    double minf;
    F_void void_stuff(&settings, &pre_calc_data, &mle_data);
    std::vector<double> param (1, SMALLTOLERANCE);
    nlopt::result result;
    double second_der, error;
    size_t iterations;

    if(! same){
      if(do_haploid_model){
        opt.set_max_objective(objective_func_F_haploid, &void_stuff);
        result = opt.optimize(param, minf);
        second_der = objective_func_F_second_deriv_haploid(param[0], pre_calc_data, mle_data);
      } else {
        opt.set_max_objective(objective_func_F, &void_stuff);
        result = opt.optimize(param, minf);
        second_der = objective_func_F_second_deriv(param[0], pre_calc_data, mle_data);
      }
      error = 1.96/std::sqrt(-second_der);
      iterations=void_stuff.iteration;

      // assign params
      last_minf = minf;
      last_param = param;
      last_result = result;
      last_error = error;
      last_iterations = iterations;
    } else {
      minf = last_minf;
      param = last_param;
      result = last_result;
      error = last_error;
      iterations = last_iterations;
    }

    // // keep this
    // double minf;
    // F_void void_stuff(settings, &pre_calc_data, &mle_data);
    // std::vector<double> param (1, 0.5);
    // nlopt::result result;
    // double second_der, error;
    // size_t iterations;
    // opt.set_max_objective(objective_func_F, &void_stuff);
    // result = opt.optimize(param, minf);
    // second_der = objective_func_F_second_deriv(param[0], pre_calc_data, mle_data);
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
    // double first_dinucl_genotype = std::exp(loglikelihood_after_mle(param, &pre_calc_data[curr_idx]));
    // double sum_exp = first_dinucl_genotype;
    // for (auto & val : pre_calc_data[curr_idx].remaining_dinucl_genotypes){
    //   sum_exp += std::exp(val);
    // }
    // end

    f << "contig: " << settings.chrom
      << " Center_pos: " << mle_data.curr_pos
      << " N_CpGs: " << mle_data.n_cpgs
      << " Depth: " << mle_data.total_depth
      << " Distance: " << mle_data.max_pos - mle_data.min_pos
      << " ll: " << minf
      << " f: " << param[0]
      << " f(95%conf): " << param[0]-error  << "," <<  param[0]+error
      << " iterations: " << iterations
      << " optim_return_code: " << result
      << " Incl_pos: " << mle_data.positions[0];
      for (auto i=mle_data.positions.begin()+1; i!=mle_data.positions.end(); i++){
        f << ","<< *i;
      }

      // // printing probability of dinucleotide genotypes for the site in the center.
      // f << " dinucl_genos: " << first_dinucl_genotype/sum_exp ;
      // for (auto & val : pre_calc_data[curr_idx].remaining_dinucl_genotypes){
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
                 std::vector<pre_calc_per_site> & pre_calc_data,
                 std::vector<std::pair<size_t, size_t>> & bed_coord) {
  std::string filename = settings.outbase + ".BED.F";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  // nlopt::opt opt(nlopt::LN_BOBYQA, 1);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  std::vector<double> lb(1, SMALLTOLERANCE), up(1, 1-SMALLTOLERANCE);
  opt.set_maxeval(1000);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(up);
  opt.set_xtol_abs(0);
  opt.set_ftol_abs(1e-15);
  opt.set_xtol_rel(0);
  opt.set_ftol_rel(0);
  for (const auto & bed : bed_coord){
    std::vector<size_t> sites_to_include;
    for (size_t curr_idx = 0; curr_idx < pre_calc_data.size(); curr_idx++) {

      if (pre_calc_data[curr_idx].position >= bed.first && pre_calc_data[curr_idx].position < bed.second){

        sites_to_include.push_back(curr_idx);
      }
      if(pre_calc_data[curr_idx].position >= bed.second){
        break;
      }

    }
    if( sites_to_include.size()==0){
      f << settings.chrom << ":" << bed.first<<"-"<<bed.second << " " << "NO_SITES_IN_THE_REGION NOT_USED" << '\n';
      continue;
    }

    // it is stupid to run over data again, but too tired to make new structs for the BED setup.
    // FIXME: make a struct that does not need to be initialized, cause then it can be used in the for loop above.
    per_mle_run mle_data(pre_calc_data[sites_to_include[0]], sites_to_include[0]);
    for (auto it=sites_to_include.begin()+1; it!=sites_to_include.end();it++){
      update_mle_data(mle_data, pre_calc_data[*it], *it);
    }

    if(mle_data.total_depth == 0){
      f << settings.chrom << ":" << bed.first<<"-"<<bed.second << " " << "NO_DATA_IN_THE_REGION NOT_USED"  << '\n';
      continue;
    }

    double minf;
    F_void void_stuff(&settings, &pre_calc_data, &mle_data);
    // std::vector<double> param (1, 0.5);
    std::vector<double> param (1, SMALLTOLERANCE);
    nlopt::result result;
    double second_der, error;
    size_t iterations;

    // // temp
    // std::vector<double> dummy;
    // for (double i=0; i<=1; i+=0.01){
    //   std::vector<double> x(i);
    //   std::cout << i << " " << std::setprecision(10) << objective_func_F(param, dummy, &void_stuff) << std::endl;
    // }
    // loglikelihood_after_mle(param[0],
    // // temp end

    if(do_haploid_model){
      opt.set_max_objective(objective_func_F_haploid, &void_stuff);
      result = opt.optimize(param, minf);
      second_der = objective_func_F_second_deriv_haploid(param[0], pre_calc_data, mle_data);
    } else {
      opt.set_max_objective(objective_func_F, &void_stuff);
      result = opt.optimize(param, minf);
      second_der = objective_func_F_second_deriv(param[0], pre_calc_data, mle_data);
    }
    error = 1.96/std::sqrt(-second_der);
    iterations=void_stuff.iteration;

    f << settings.chrom << ":" << bed.first<<"-"<<bed.second
      << " N_CpGs: " << mle_data.n_cpgs
      << " Depth: " << mle_data.total_depth
      << " Distance: " << mle_data.max_pos - mle_data.min_pos
      << " ll: " << minf
      << " f: " << param[0]
      << " f(95%conf): " << param[0]-error  << "," <<  param[0]+error
      << " iterations: " << iterations
      << " optim_return_code: " << result
      << " Incl_pos: " << mle_data.positions[0];
      for (auto i=mle_data.positions.begin()+1; i!=mle_data.positions.end(); i++){
        f << ","<< *i;
      }
      f << std::endl;

  }
  f.close();
}

double base_condition_one_genotype(const size_t &strand, const double &deam,
                                   const double &seqerror, const size_t &base,
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

double calc_prob_obs_base(const double & seqerror, const size_t & base_compos, const double & deamin_rate){
  double res;
  if(base_compos == 0){  // CC :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
    res = ((1-deamin_rate) * (1-seqerror)) + (deamin_rate * seqerror/3) + (seqerror * seqerror/3);
  }else if(base_compos == 1){  // CT :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
    res = (deamin_rate * (1-seqerror)) + ((1-deamin_rate) * seqerror/3) + (seqerror * seqerror/3);
  } else { // mutation
    res = ((1-deamin_rate) * seqerror/3) + (deamin_rate * seqerror/3) + (seqerror * seqerror/3);
  }
  return res;
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
    std::cerr << "\t -> ERROR: provide a string of seven prior e.g. -h 1,0,0,0,0,0,0" << '\n';
    exit(EXIT_FAILURE);
  }

  // for (const auto & val : LOG_PRIORS){
  //   std::cerr << val << '\n';
  // }

  for (size_t i=0; i<res.size(); i++){
    LOG_PRIORS[i] = std::log(res[i]/sum);
  }

  // for (const auto & val : LOG_PRIORS){
  //   std::cerr << val << '\n';
  // }

}


void parse_bed_file(general_settings & settings, std::vector<std::pair<size_t, size_t>> & res){
  std::string filename = settings.bed_f;
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    std::string row;
    std::string chrom;
    size_t start,end ;
    std::stringstream ss;
    while(getline(f, row)){
      if(row.empty()){
        std::cerr << "BED: "<< settings.bed_f << " contains empty rows. EXITING" << '\n';
        exit(EXIT_FAILURE);
      }
      ss.str(row);

      ss >> chrom >> start >> end;
      if(start>=end){
        std::cerr << "BED: "<< settings.bed_f << " is not normal: chr: " << chrom << " Start: " << start << " END: " << end << ". EXITING" << '\n';
        exit(EXIT_FAILURE);
      }
      if(chrom==settings.chrom){
        res.push_back(std::make_pair (start,end));
      }
      ss.clear();
    }
  }
  f.close();
}

void print_d(pre_calc_per_site & d){
  std::cout << d.depth << " " << d.position <<  '\n';
  for (const auto & val : d.remaining_dinucl_genotypes){
    std::cout << val <<  '\n';
  }
  for (int i=0; i<d.depth;i++){
    std::cout << i << " " << d.maperrors[i] << " " << d.pre_noM[i] << " " << d.pre_M[i] << '\n';
  }
  std::cout << std::flush;
}

void print_help(){
  std::cerr << "\nDamMet (" << DAMMET_VERSION << ") is a software aimed to estimate methylation maps using HTS sequencing "
    "data underlying ancient samples. The implemented model follows a two-steps procedure. "
    "The first step obtains a Maximum Likelihood Estimate (MLE) of position-specific deamination "
    "rates at both methylated and unmethylated cytosine residues. The second step makes use "
    "of these estimates to recover a MLE of local methylation levels in a user-defined window size." << std::endl;
  std::cerr << "Three args are required:\n\t->-b (bam)\n\t->-r (reference "
    "fasta)\n\t->-c (chromosome of interest)" << '\n';
  std::cerr << "OPTIONS:" << std::endl;
  std::cerr << "\t-> BED file (-B): " << std::endl;
  std::cerr << "\t-> minmapQ (-q): " << std::endl;
  std::cerr << "\t-> minbaseQ (-Q): " << std::endl;
  std::cerr << "\t-> MinReadLength (-L): " << std::endl;
  std::cerr << "\t-> MinReadLength_Deamrates (-l): " << std::endl;
  std::cerr << "\t-> Max_Pos_From_End (-P): " << std::endl;
  std::cerr << "\t-> Expected fraction of methylated CpGs (-M): " << std::endl;
  std::cerr << "\t-> Outbase (-O): " << std::endl;
  std::cerr << "\t-> readFlags (-F): " << std::endl;
  std::cerr << "\t-> Number of cycles (-C) (Only used if no RG file is NOT provided): " << std::endl;
  std::cerr << "\t-> Using Precalc deamination rates from (-D): " << std::endl;
  std::cerr << "\t-> Using readgroups from file (-R): " << std::endl;
  std::cerr << "\t-> Exclude sites (1-based) (-E): " << std::endl;
  std::cerr << "\t-> Exclude BED (-e): " << std::endl;
  std::cerr << "\t-> WindowSize (-W): " << std::endl;
  std::cerr << "\t-> Max CpGs per Window (-N): " << std::endl;
  exit(EXIT_SUCCESS);
}


void est_dam_only(general_settings & settings) {
  std::string stream_filename (settings.outbase+".args");
  settings.args_stream.open(stream_filename.c_str());
  checkfilehandle(settings.args_stream, stream_filename);
  settings.args_stream << settings.all_options;

  std::vector<std::pair<size_t, size_t>> bed_coord;
  // parsing bedfile if provided
  if(!settings.bed_f.empty()){
    parse_bed_file(settings, bed_coord);
    settings.buffer += "\t-> Analyzing: " + std::to_string(bed_coord.size()) + " BED regions on chrom: " + settings.chrom + '\n';
    print_log(settings);

    if(bed_coord.size()==0){
      std::cerr << '\n' << "EXITING. NOT BED COORDS ON CHROM" << '\n';
      exit(EXIT_SUCCESS);
    }
  }

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

  // Load the index
  hts_idx_t *idx = sam_index_load(in, settings.bam_fn.c_str());

  if (idx == NULL) {
    std::cerr
        << "Unable to open BAM/SAM index. Make sure alignments are indexed."
        << '\n';
    exit(EXIT_FAILURE);
  }

  // load read groups
  bool rg_split;
  std::vector<std::string> rgs;
  std::vector<size_t> cycles;
  if((!settings.readgroups_f.empty()) && check_file_exists(settings.readgroups_f)){
    rg_split = true;
    std::ifstream f (settings.readgroups_f.c_str());
    checkfilehandle(f, settings.readgroups_f);
    std::string row, rg, cycle_string;
    size_t cycle;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> rg >> cycle_string;
      if(cycle_string.empty()){
        cycle=std::numeric_limits<size_t>::max();
        settings.buffer += "\t-> RG: " + rg + '\n';
        print_log(settings);
      } else {
        cycle=std::stoi(cycle_string);
        settings.buffer += "\t-> RG: " + rg + ". Sequencing cycles: " + cycle_string + '\n';
        print_log(settings);
      }
      rgs.push_back(rg);
      cycles.push_back(cycle);
      ss.clear();
      rg.clear();
      cycle_string.clear();
    }
    f.close();
  } else {
    rg_split = false;
    rgs.push_back(ALL_RG);
    settings.buffer += "\t-> Merging all reads into a single read group named: " + ALL_RG;
    if(settings.cycles!=std::numeric_limits<size_t>::max()){
      cycles.push_back(settings.cycles);
      settings.buffer += ". Sequencing cycles: " + std::to_string(settings.cycles);
    } else {
      cycles.push_back(std::numeric_limits<size_t>::max());
    }
    settings.buffer += '\n';
    print_log(settings);
  }

  // go to the correct chromosome
  // https://github.com/gatoravi/bam-parser-tutorial/blob/master/parse_bam.cc
  hts_itr_t *iter = sam_itr_querys(idx, header, settings.chrom.c_str());

  // load ref of that chromosome
  faidx_t *fai = ref_init(settings.reference_fn);

  char *ref = fetch_chrom(fai, settings.chrom);

  size_t seq_len = faidx_seq_len(fai, settings.chrom.c_str());

  if (!settings.exclude_bed_fn.empty()) {
    settings.buffer += "\t-> Masking genomic BED regions (-e) " + settings.exclude_bed_fn + '\n';
    print_log(settings);
    filter_ref_bed(settings.chrom, settings.exclude_bed_fn, ref);
  }

  if (!settings.exclude_sites_fn.empty()) {
    settings.buffer += "\t-> Masking genomic sites sites (-E) " + settings.exclude_sites_fn + '\n';
    print_log(settings);
    filter_ref_sites(settings.chrom, settings.exclude_sites_fn, ref);
  }

  std::vector<std::vector<double>> tmf;
  tmf.resize(rgs.size());
  std::vector<std::vector<int>> tm;
  tm.resize(rgs.size());
  for (size_t i=0; i<rgs.size(); i++){
    tm[i] = init_tallymat<int>(settings.max_pos_to_end);
  }

  my_cov_rg cov_rg(rgs.size());

  const keeplist_map cpg_map = get_cpg_chrom_pos(ref, seq_len);
  // const std::vector<int> cpg_bool = get_cpg_chrom_bool(ref, seq_len);
  const v_un_ch cpg_bool = get_cpg_chrom(ref, seq_len);
  settings.buffer += "\t-> " + std::to_string(cpg_map.size()) + " CpG's in chrom: " + settings.chrom + '\n';
  print_log(settings);

  // do not included CnonCpGs if deamrates are already provided.
  std::vector<size_t> exclude_CnonCpGs(rgs.size(), 0);

  std::vector<uni_ptr_obs> cpg_data, nocpg_data;
  cpg_data.resize(rgs.size());
  nocpg_data.resize(rgs.size());
  // std::vector<per_site> data;
  // data.resize(cpg_map.size());
  //size_t nocpg_to_include = 1e6;
  // per_site_nocpg nocpg_data(nocpg_to_include);

  bam1_t *rd = bam_init1();
  alignment_data d;
  int reads = 1e6;
  int counter = 0;

  size_t trashed=0, reads_skipped=0;
  size_t startpos, endpos;
  bool skipread=true;
  time_t start_time_load_data,end_time_load_data;

  uint8_t *rgptr;
  std::string rgname;
  size_t rgname_idx;
  std::string rname;
  time(&start_time_load_data);
  while (sam_itr_next(in, iter, rd) >= 0) {
    skipread = true;
    if (counter % reads == 0 && reads) {
      std::cerr << "\t-> " << counter << " reads processed and " << trashed << " discarded. " << '\r';
    }
    counter++;

    // check that a cpg is present
    // this is 50% faster than analyzing every single read.
    startpos = rd->core.pos;
    endpos = bam_endpos(rd);
    for (size_t i=startpos; i<=endpos; i++){
      if(cpg_bool[i]==1){
        skipread=false;
        break;
      }
    }

    if(skipread){
      reads_skipped++;
      continue;
    }

    if (rd->core.l_qseq < settings.minreadlength ||
        rd->core.qual < settings.minmapQ ||
        (rd->core.flag & settings.flags_off) != 0 ||
        rd->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
      trashed++;
      continue;
    }

    if (rg_split) {
      rgptr = bam_aux_get(rd, "RG");
      if (rgptr == NULL) {
        // in the case norg exists. dont know if i should trash or put in
        // UNKWOWN
        rname = std::string(bam_get_qname(rd));
        std::cerr << "[No_found_RG] :: " << rname << '\n';
        trashed++;
        continue;
      } else {
        // moves past the Z by +1. all
        // colons are gone already
        rgname = std::string((const char *)(rgptr + 1));
        auto match = std::find(rgs.begin(), rgs.end(), rgname);
        if (match != rgs.end()) {
          rgname_idx = distance(rgs.begin(), match);
        } else { // in the case it cannot find the rg. The read will be trashed
          trashed++;
          continue;
        }
      }
   } else {
     rgname=ALL_RG;
     rgname_idx=ALL_DEAMMETH_RG_IDX;
   }
   d = align_read(rd, ref);
   add_aligned_data(settings, d,
                    cpg_data[rgname_idx],
                    nocpg_data[rgname_idx],
                    tm[rgname_idx],
                    cov_rg, rgname_idx,
                    cycles[rgname_idx],
                    exclude_CnonCpGs[rgname_idx]);
  }
  time(&end_time_load_data);

  settings.buffer += "\t-> Processed: " + std::to_string(counter)  +
    ". Reads filtered: " + std::to_string(trashed) +
    ". Reads skipped (nocpg overlap): " + std::to_string(reads_skipped) +
    ". Loaded in " + std::to_string(difftime(end_time_load_data, start_time_load_data)) + " seconds." + '\n';


  print_log(settings);
  for (size_t i=0; i<rgs.size(); i++){
    settings.buffer += "\t-> Total_Observations RG: "+ rgs[i] +
      " CpG/cnonCpG: " + std::to_string(cov_rg.cpg[i]) + " " + std::to_string(cov_rg.nocpg[i]) +
      " CpGcov: " + std::to_string((double)cov_rg.cpg[i]/(double)cpg_map.size()) + '\n';
    print_log(settings);
  }

  // dumping counts to file
  for (size_t i=0; i<rgs.size(); i++){
    settings.buffer += "\t-> Dumping count file: " + settings.outbase + "." + rgs[i] +  ".tallycounts" + '\n';
    print_log(settings);
    dump_count_file(settings, tm[i], rgs[i]);
    tmf[i] = read_count_file(settings, rgs[i]);
  }

  std::vector<std::vector<double>> param_deam;
  param_deam.resize(rgs.size());

  settings.args_stream << std::flush;
  for (size_t i=0; i<rgs.size(); i++){
    // if((!settings.deamrates_filename.empty()) && check_file_exists(settings.deamrates_filename)){
    //   settings.buffer += "\t-> Loading deamination rates from " + settings.deamrates_filename + '\n';
    //   print_log(settings);
    //   std::cerr << "\t-> Make sure that the file contains the same number of pos to include. Deammeth does not check that" << '\n';
    //   param_deam[i] = load_deamrates_f(settings);
    // } else if (check_file_exists(settings.outbase+"."+rgs[i]+".deamrates")){
    //   settings.buffer += "\t-> Loading deamination rates from " + settings.outbase+"."+rgs[i]+".deamrates" + '\n';
    //   print_log(settings);
    //   param_deam[i] = load_deamrates(settings, rgs[i]);
    // }else {
      settings.buffer += "\t-> Starting Optim of deamination rates. RG: " + rgs[i] + '\n';
      print_log(settings);
#if 1
      print_data(cpg_data[i], nocpg_data[i]);
#endif
      run_deamrates_optim(settings, tmf[i], cpg_data[i], nocpg_data[i], rgs[i]);
      settings.buffer += "\t-> Dumping deamination parameters to " + settings.outbase+"."+rgs[i]+".deamrates" + '\n';
      print_log(settings);
      param_deam[i] = load_deamrates(settings, rgs[i]);
      // }
  }
  std::cerr << "\t-> Cleaning up." << '\n';
  settings.args_stream << std::flush;
  free(ref);
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(rd);
  bam_hdr_destroy(header);
  sam_close(in);
}

int est_dam_and_F(general_settings & settings){
  return 1;
}

int main(int argc, char *argv[]) {
  time_t start_time, end_time;
  time(&start_time);
  general_settings settings;
  args_parser(argc, argv, settings);

  if (settings.bam_fn.empty() || settings.reference_fn.empty() ||
      settings.chrom.empty()) {
    print_help();
  }
  if(settings.max_cpgs==std::numeric_limits<size_t>::max() && settings.windowsize==std::numeric_limits<size_t>::max() && settings.bed_f.empty()){
    std::cerr << "Must specify either -N (max CpGs per window) AND/OR -W (max windowsize) OR -B bedfile" << '\n';
    std::cerr << "EXITING...." << '\n';
    exit(EXIT_FAILURE);
  }


  // last option should be to get F from est dam only
  if(0){
    int a=est_dam_and_F(settings);
  } else {
    est_dam_only(settings);
  }
  time(&end_time);
  double time_gone = difftime(end_time, start_time);
  size_t minutes = time_gone / 60.0;
  size_t seconds = (int)time_gone % 60;
  settings.buffer += "\t-> Done in " + std::to_string(minutes)+ ":" + std::to_string(seconds) + " M:S." + '\n';
  print_log(settings);
}
