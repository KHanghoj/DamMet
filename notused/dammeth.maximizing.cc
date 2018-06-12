#include "dammeth.h"
#include "nlopt.hpp"

int break_MCMC=0;

template <typename T>
void checkfilehandle(T &fh, std::string filename){
  if (! fh.is_open()){
    std::cerr << "Couldnt open file: " << filename << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
}

bool check_file_exists(std::string filename){
  std::ifstream f(filename.c_str(), std::ios::in);
  return f.good();
}

double phred_to_double(const int & phred){
  return std::pow(10, -(phred/10.0));
}

// readpos, prime, strand, refdinucl1, refdinucl2, sampledinucl1,sampledinucl2
using tally_mat = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<size_t>>>>>>>;
using tally_mat_freq = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<double>>>>>>>;

tally_mat init_tallymat(const size_t readpos){
  tally_mat thematrix;
  thematrix.resize(readpos+2);
  for (size_t p=0; p<(readpos+2); p++ ){
    thematrix[p].resize(PRIMES);

    for (size_t pr=0; pr<PRIMES; pr++){
      thematrix[p][pr].resize(STRANDS);

      for (size_t strand=0; strand < STRANDS; strand++){
        thematrix[p][pr][strand].resize(NUCLEOTIDES);

        for (size_t r1=0; r1<NUCLEOTIDES; r1++){
          thematrix[p][pr][strand][r1].resize(NUCLEOTIDES);

          for (size_t r2=0; r2<NUCLEOTIDES; r2++){
            thematrix[p][pr][strand][r1][r2].resize(NUCLEOTIDES);

            for (size_t s1=0; s1<NUCLEOTIDES; s1++){
              thematrix[p][pr][strand][r1][r2][s1].resize(NUCLEOTIDES);
            }
          }
        }
      }
    }
  }
  return thematrix;
}

tally_mat_freq init_tallymat_freq(const size_t readpos){
  tally_mat_freq thematrix;
  thematrix.resize(readpos+2);
  for (size_t p=0; p<(readpos+2); p++ ){
    thematrix[p].resize(PRIMES);

    for (size_t pr=0; pr<PRIMES; pr++){
      thematrix[p][pr].resize(STRANDS);

      for (size_t strand=0; strand < STRANDS; strand++){
        thematrix[p][pr][strand].resize(NUCLEOTIDES);

        for (size_t r1=0; r1<NUCLEOTIDES; r1++){
          thematrix[p][pr][strand][r1].resize(NUCLEOTIDES);

          for (size_t r2=0; r2<NUCLEOTIDES; r2++){
            thematrix[p][pr][strand][r1][r2].resize(NUCLEOTIDES);

            for (size_t s1=0; s1<NUCLEOTIDES; s1++){
              thematrix[p][pr][strand][r1][r2][s1].resize(NUCLEOTIDES);
            }
          }
        }
      }
    }
  }
  return thematrix;
}


// code for MCMC estimates of deamination for unmeth CpG and meth CpG
// doing it for plus strand only
int D_LOWER = 0;
int D_UPPER = 1;

// tmf: // readpos, prime, strand, refdinucl1, refdinucl2, sampledinucl1,sampledinucl2
// meth/unmeth, readpos, prime
using mcmc_meth_mat = std::vector<std::vector<std::vector< double > > >;

mcmc_meth_mat get_vectors(const general_settings & settings){
  mcmc_meth_mat thematrix;
  // 0 is meth; 1 is unmeth.
  thematrix.resize(2);
  for (size_t i=0; i<METHSTATES;i++){
    thematrix[i].resize(settings.max_pos_to_end+2);
    for (size_t p=0; p<settings.max_pos_to_end+2;p++){
      thematrix[i][p].resize(PRIMES);
    }
  }
  return thematrix;
}

void merge_vectors(const general_settings & settings, mcmc_meth_mat & mat_to, mcmc_meth_mat & mat_from){
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<settings.max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
        // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
        mat_to[i][p][pr] += mat_from[i][p][pr];
      }
    }
  }
}

void normalize_vectors(const general_settings & settings, mcmc_meth_mat & final_mat, size_t & mat_aggregated){
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<settings.max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
        // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
        final_mat[i][p][pr] /= (double)mat_aggregated;
      }
    }
  }
}

mcmc_meth_mat init_mcmc_meth_mat(const general_settings & settings, const tally_mat_freq & tmf){
  mcmc_meth_mat thematrix = get_vectors(settings);
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<settings.max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
        // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
        // double val = dist_normal(gen);
        // while (val > D_UPPER && val < D_LOWER){
        //   val = dist_normal(gen);
        // }
        // thematrix[i][p][pr] = val;

        // using the count matrix
        if(i==0){ // methylated cpg
          thematrix[i][p][pr] = tmf[p][pr][0][1][2][3][2];
        } else {
          thematrix[i][p][pr] = tmf[p][pr][0][1][1][3][1];
        }

        // for simple 10 read simulations put all deam values to 0:
        // thematrix[i][p][pr] = 0;
      }
    }
  }
  return thematrix;
}

mcmc_meth_mat init_zero_mcmc_meth_mat(const general_settings & settings){
  mcmc_meth_mat thematrix = get_vectors(settings);
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<settings.max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
        // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
        thematrix[i][p][pr] = 0;
      }
    }
  }
  return thematrix;
}


const keeplist_map get_cpg_genomic_pos(const char * ref, const size_t & seq_len){
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

void dump_count_file(general_settings & settings, tally_mat tm){
  std::string filename = settings.outbase + ".tallycounts";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  if (f.is_open()) {
    size_t readpos = settings.max_pos_to_end;
    for (size_t p = 0; p < (readpos + 2); p++) {
      for (size_t pr = 0; pr < PRIMES; pr++) {
        for (size_t strand = 0; strand < STRANDS; strand++) {
          for (size_t r1 = 0; r1 < NUCLEOTIDES; r1++) {
            for (size_t r2 = 0; r2 < NUCLEOTIDES; r2++) {
              for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
                for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                  if (tm[p][pr][strand][r1][r2][s1][s2]) {
                    f << p << " " << pr << " " << strand << " " << r1 << " "
                      << r2 << " " << s1 << " " << s2 << " "
                      << tm[p][pr][strand][r1][r2][s1][s2] << '\n';
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  f.flush();
  f.close();
}

tally_mat_freq read_count_file(general_settings & settings){
  tally_mat temp = init_tallymat(settings.max_pos_to_end);
  tally_mat_freq res = init_tallymat_freq(settings.max_pos_to_end);
  std::string filename = settings.outbase + ".tallycounts";
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    std::string row;
    size_t p, pr, strand, r1, r2, s1, s2;
    double val;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> p >> pr >> strand >> r1 >> r2 >> s1 >> s2 >> val;
      temp[p][pr][strand][r1][r2][s1][s2] = val;
      ss.clear();
    }
  }
  f.close();

  size_t readpos = settings.max_pos_to_end;
  for (size_t p = 0; p < (readpos + 2); p++) {
    for (size_t pr = 0; pr < PRIMES; pr++) {
      for (size_t strand = 0; strand < STRANDS; strand++) {
        for (size_t r1 = 0; r1 < NUCLEOTIDES; r1++) {
          for (size_t r2 = 0; r2 < NUCLEOTIDES; r2++) {
            double totalsum = 0;
            for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
              for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                if (temp[p][pr][strand][r1][r2][s1][s2]) {
                  totalsum += temp[p][pr][strand][r1][r2][s1][s2];
                }
              }
            }
            for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
              for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                if (temp[p][pr][strand][r1][r2][s1][s2]) {
                  res[p][pr][strand][r1][r2][s1][s2] = temp[p][pr][strand][r1][r2][s1][s2] / totalsum;

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

// random device class instance, source of 'true' randomness for initializing random seed
std::random_device rd;
// Mersenne twister PRNG, initialized with seed from previous random device instance
std::mt19937 gen(rd());

void print_mcmc_meth_mat(const general_settings & settings, std::ofstream & f , mcmc_meth_mat & d, const size_t & chain, long double & ll ){
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<settings.max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
        // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
#ifdef skipinternalpos
        if (p==settings.max_pos_to_end+1){
          f << chain << " " << i << " " << p << " " << pr << " " << d[i][p-1][pr] << " " << ll<< '\n';
        } else {
          f << chain << " " << i << " " << p << " " << pr << " " << d[i][p][pr] << " " << ll<< '\n';
        }
#else
        f << chain << " " << i << " " << p << " " << pr << " " << d[i][p][pr] << " " << ll<< '\n';
#endif

      }
    }
  }
  f << std::flush ;

}

void get_new_random_values(const general_settings &settings,
                           const mcmc_meth_mat &d_prev, mcmc_meth_mat &d_next,
                           std::normal_distribution<double> &dist_normal) {
  for (size_t i=0; i<METHSTATES;i++){  // 0 is meth; 1 is unmeth.
    for (size_t p=0; p<(settings.max_pos_to_end+2);p++){
      for (size_t pr=0; pr<PRIMES; pr++){
        if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
        // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
        // std::normal_distribution<double> dist_normal (d_prev[i][p][pr], DIST_DEVIATION);
        double val = d_prev[i][p][pr] + dist_normal(gen);
        if (val > D_LOWER && val < D_UPPER){
          d_next[i][p][pr] = val;
        } else {
          d_next[i][p][pr] = d_prev[i][p][pr];
        }
      }
    }
  }
}

// Printing dinucl_data
void print_dinucl(const general_settings &settings,
                  const std::vector<per_site> &data,
                  const tally_mat_freq &tmf) {

  std::string filename = settings.outbase + ".dinucls";
  std::ofstream f(filename.c_str());
  checkfilehandle(f, filename);
  if (f.is_open()) {
    for (const auto &site : data) {
      dinucl_pair_of_pair *dinucl_max =
        &ALL_DINUCL_PAIR_OF_PAIRS.at(site.dinucl_max_idx);
      dinucl_pair_of_pair *t = &ALL_DINUCL_PAIR_OF_PAIRS.at(IDX_CpG_CpG);
      f << site.position << " " << site.depth << " " << t->first.first
        << t->first.second << " " << t->second.first << t->second.second << " "
        << site.dinucl_genotype_likelihoods.at(IDX_CpG_CpG) << " "
        << site.dinucl_max_idx << " "
        << site.dinucl_genotype_likelihoods[site.dinucl_max_idx] << " "
        << dinucl_max->first.first << dinucl_max->first.second << " "
        << dinucl_max->second.first << dinucl_max->second.second << " ::: add simple count statistic 0,1,2";
      f << '\n';
    }
  }
  f.close();
}



void get_bases(const alignment_data & d, const int & b1, const int & b2, int & r_base1, int & r_base2, int & s_base1, int & s_base2){
  r_base1 = d.t_ref[b1];
  r_base2 = d.t_ref[b2];
  s_base1 = d.t_seq[b1];
  s_base2 = d.t_seq[b2];
}

void get_pos_and_prime(const general_settings & settings, const size_t & dist_3p, const size_t & dist_5p, size_t & pos, size_t & prime){
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
  // if (dist_5p <= dist_3p){
  //   if(dist_5p>settings.max_pos_to_end){
  //     pos = settings.max_pos_to_end + 1;
  //     prime = 0; // 5p+1 by default
  //   } else {
  //     pos = dist_5p;
  //     prime = 0; // 5p+1 by default
  //   }
  // } else {
  //   if(dist_3p>settings.max_pos_to_end){
  //     pos = settings.max_pos_to_end + 1;
  //     prime = 1; // 5p+1 by default
  //   } else {
  //     pos = dist_3p;
  //     prime = 1; // 5p+1 by default
  //   }
  // }
}

bool nucleotide_missing(const int & r_base1, const int & r_base2, const int & s_base1, const int & s_base2){
  return (r_base1==4 || r_base2==4 || s_base1==4 || s_base2==4);
}

alignment_data align_read(bam1_t * rd, char * ref){

  alignment_data d ;
  d.n_nucleotides = rd->core.l_qseq;
  d.strand = bam_is_rev(rd);
  d.mapQ = rd->core.qual;
  if(d.mapQ>37){
    d.mapQ = 37;
  }

  uint8_t *quals = bam_get_qual(rd);
  uint8_t *seq = bam_get_seq(rd);
  int nCig = rd->core.n_cigar;
  uint32_t *cigs = bam_get_cigar(rd);
  int seq_pos = 0;              // position within sequence
  int wpos = rd->core.pos; // this value is the current position
                                // assocatied with the positions at seq_pos
  int hasInfo = 0;
  int c;

  // loop through read by looping through the cigar string
  for (int i = 0; i < nCig; i++) {
    int opCode = cigs[i] & BAM_CIGAR_MASK;  // what to do
    int opLen = cigs[i] >> BAM_CIGAR_SHIFT; // length of what to do
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
          d.t_ref.push_back(4);

          d.t_positions.push_back(0);

          seq_pos++; // <- important, must be after macro
        }
        wpos++;
      } else { // this is the deletion part
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

void tally_aligned_data(const general_settings & settings, const alignment_data & d, tally_mat & tm, char * ref_tally){
  int b1, b2, r_base1, r_base2, s_base1, s_base2;
  size_t dist_5p, dist_3p;
  size_t pos, prime;
  for (size_t i = 0; i < d.t_seq.size()-1; i++){
    if(d.strand){ // negative strand
      b1 = d.t_seq.size()-i-2;
      b2 = d.t_seq.size()-i-1; // closest to 5 prime
      get_bases(d, b1, b2, r_base1, r_base2, s_base1, s_base2);
      dist_5p = d.t_isop[b2];
      dist_3p = d.t_posi[b2];

    } else { // positive strand
      b1 = i; // closest to 5 prime
      b2 = i+1;
      get_bases(d, b1, b2, r_base1, r_base2, s_base1, s_base2);
      dist_5p = d.t_posi[b1];
      dist_3p = d.t_isop[b1];
    }

    if( nucleotide_missing(r_base1, r_base2, s_base1, s_base2) ||  d.t_qs[b1] < settings.minbaseQ || d.t_qs[b2] < settings.minbaseQ ){
      continue;
    }

    if(refToInt[(int)ref_tally[d.t_positions[b1]]]==4 || refToInt[(int)ref_tally[d.t_positions[b2]]]==4){
      continue;
    }

    get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

    // discard 3prime data if read == number of cycles. as we do not know the actual end of the read.
    if (prime==1 && d.n_nucleotides==settings.cycles){
      continue;
    }

    tm[pos][prime][d.strand][r_base1][r_base2][s_base1][s_base2]++;
  }
}

bool no_ref_cpg(const int & r_base1, const int & r_base2){
  return(r_base1!=1 || r_base2!=2);
}

bool cpg(const int & base1, const int & base2){
  return(base1==1 && base2==2);
}

void add_aligned_data(const general_settings &settings, const alignment_data &d,
                      std::vector<per_site> &data, const keeplist_map &cpg_map,
                      per_site_nocpg &nocpg_data, const char *ref_tally) {
  per_site *p_per_site;
  int b1, b2, r_base1, r_base2, s_base1, s_base2;
  size_t dist_5p, dist_3p;
  size_t pos, prime;
  size_t data_idx=0;

  // MCMC:
  double seqerror = 0;
  double maperror = phred_to_double(d.mapQ);
  int base_composition; // 0,1,2

  for (size_t i = 0; i < d.t_seq.size()-1; i++){
    if(d.strand){ // negative strand
      b1 = d.t_seq.size()-i-2;
      b2 = d.t_seq.size()-i-1; // closest to 5 prime
      get_bases(d, b1, b2, r_base1, r_base2, s_base1, s_base2);
      if(nucleotide_missing(r_base1, r_base2, s_base1, s_base2) ||  d.t_qs[b1] < settings.minbaseQ || d.t_qs[b2] < settings.minbaseQ ){
        continue;
      }

      dist_5p = d.t_isop[b2];
      dist_3p = d.t_posi[b2];
      seqerror = phred_to_double(d.t_qs[b2]);
      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

#ifdef skipinternalallanal
      // // NOTE: this removes all internal sites. Not sure it is a good idea in general
      if (prime==0 && pos == (settings.max_pos_to_end+1)){
         continue;
      }
#endif
      // discard 3prime data, if read == number of cycles.
      if (prime==1 && d.n_nucleotides==settings.cycles){
        continue;
      }


      if (s_base2 == 2){
        base_composition = 0;  // no error
      } else if (s_base2 == 0){
        base_composition = 1;  // deamin or error
      } else {
        base_composition = 2; // error (could also be deamin then error)
      }

      if(no_ref_cpg(r_base1, r_base2)){
        if(refToInt[(int)ref_tally[d.t_positions[b2]]]==2 && nocpg_data.depth < nocpg_data.max_depth){
          nocpg_data.pos_to_end.push_back(pos);
          nocpg_data.prime.push_back(prime);
          nocpg_data.base_compos.push_back(base_composition);
          nocpg_data.seqerrors.push_back(seqerror);
          nocpg_data.seqerrors.push_back(maperror);
          nocpg_data.depth++;
        }
        continue;
      }

    } else { // positive strand
      b1 = i; // closest to 5 prime
      b2 = i+1;
      get_bases(d, b1, b2, r_base1, r_base2, s_base1, s_base2);
      if(nucleotide_missing(r_base1, r_base2, s_base1, s_base2) ||  d.t_qs[b1] < settings.minbaseQ || d.t_qs[b2] < settings.minbaseQ ){
        continue;
      }

      dist_5p = d.t_posi[b1];
      dist_3p = d.t_isop[b1];
      seqerror = phred_to_double(d.t_qs[b1]);

      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);
#ifdef skipinternalallanal
      // // NOTE: this removes all internal sites. Not sure it is a good idea in general
      if (prime==0 && pos == (settings.max_pos_to_end+1)){
         continue;
      }
#endif

      // discard 3prime data, if read == number of cycles.
      if (prime==1 && d.n_nucleotides==settings.cycles){
        continue;
      }


      if (s_base1 == 1){
        base_composition = 0;  // no error
      } else if (s_base1 == 3){
        base_composition = 1;  // deamin or error
      } else {
        base_composition = 2; // error (could also be deamin then error)
      }
      if(no_ref_cpg(r_base1, r_base2)){
        if(refToInt[(int)ref_tally[d.t_positions[b1]]]==1 && nocpg_data.depth < nocpg_data.max_depth){
          nocpg_data.pos_to_end.push_back(pos);
          nocpg_data.prime.push_back(prime);
          nocpg_data.base_compos.push_back(base_composition);
          nocpg_data.seqerrors.push_back(seqerror);
          nocpg_data.seqerrors.push_back(maperror);
          nocpg_data.depth++;
        }
        continue;
      }

    }
    data_idx = cpg_map.at(d.t_positions[b1]);
    p_per_site = &data[data_idx];
    p_per_site->position = d.t_positions[b1];

    p_per_site->depth++;
    p_per_site->strand.push_back(d.strand);


    p_per_site->pos_to_end.push_back(pos);
    p_per_site->prime.push_back(prime);

    p_per_site->base_compos.push_back(base_composition);
    p_per_site->seqerrors.push_back(seqerror);
    p_per_site->maperrors.push_back(maperror);

    // // this is used for dinucleotide genotype likelihoods and printing.
    p_per_site->refs.push_back(std::make_pair(r_base1, r_base2));
    p_per_site->bases.push_back(std::make_pair(s_base1, s_base2));
    p_per_site->quals.push_back(std::make_pair(d.t_qs[b1], d.t_qs[b2]));
  }
}

bool check_max_dinucl_CpG(const per_site & site, const double & min_CpG_CpG_prob){
  // IDX_CpG_CpG  // 81
  // IDX_CpG_CpA  // 60
  // IDX_CpG_TpG  // 89
  // IDX_TpG_TpG  // 133
  // IDX_CpA_CpA  // 58
  return ((site.dinucl_max_idx == IDX_CpG_CpG ||
           site.dinucl_max_idx == IDX_CpG_CpA ||
           site.dinucl_max_idx == IDX_CpG_TpG ||
           site.dinucl_max_idx == IDX_TpG_TpG ||
           site.dinucl_max_idx == IDX_CpA_CpA) &&
          (site.dinucl_genotype_likelihoods[IDX_CpG_CpG] >= min_CpG_CpG_prob)
          ? true
          : false);
}

double BASE_FREQ_FLAT_PRIOR = 0.25;

struct blabla {
  general_settings settings;
  std::vector<per_site> * data;
  per_site_nocpg * nocpg_data;
  size_t iteration;
  blabla(const general_settings & s, std::vector<per_site> * d, per_site_nocpg * nocpg_d){
    settings = s;
    data = d;
    nocpg_data = nocpg_d;
    iteration = 0;
  }
};

// long double mcmc_loglike_incl_nocpg(const general_settings &settings,
//                                     const mcmc_meth_mat &d,
//                                     const std::vector<per_site> &data,
//                                     const per_site_nocpg &nocpg_data) {

mcmc_meth_mat load_deamin_probabilities_temp(const general_settings & settings, std::string  & filename ){
  mcmc_meth_mat res = get_vectors(settings);
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    std::string row;
    size_t methstate, pos, prime;
    double rate;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> methstate >> pos >> prime >> rate;
      ss.clear();
      res[methstate][pos][prime] = rate;
    }
  }
  f.close();
  return res;

}

size_t get_param_idx(const size_t & max_pos_to_end, const size_t & meth, const size_t & pos_to_end, const size_t & prime){
  return (meth * PRIMES * ( max_pos_to_end+2 ) + pos_to_end * PRIMES  + prime);
}

double objective_func(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
  if (!grad.empty()) {
    std::cerr << "trouble" << '\n';
    // grad[0] = 0.0;
    // grad[1] = 0.5 / sqrt(x[1]);
  }
  blabla *d = static_cast<blabla*>(my_func_data);
  d->iteration++;

  double ll = 0;
  double res=0;
  double seqerror, deamin_methylated, deamin_unmethylated;
  double maperror;
  double noM, M;
#ifdef skipinternalpos
  double total_depth=0, internal_depth=0;
#endif
  for (const auto &site : (*d->data)) {
    for (size_t i=0; i < site.depth; i++){

#ifdef skipinternalpos
      // -Dskipinternalpos
      if(d->site->pos_to_end[i]==(d->settings.max_pos_to_end+1)){
        internal_depth++;
        continue;
      }
      total_depth++;
#endif
      res=0;
      seqerror = site.seqerrors[i];
      maperror = site.maperrors[i];
      deamin_methylated = x[get_param_idx(d->settings.max_pos_to_end, METHSTATE, site.pos_to_end[i], site.prime[i])];
      // deamin_methylated = d[0][site.pos_to_end[i]][site.prime[i]];
      deamin_unmethylated = x[get_param_idx(d->settings.max_pos_to_end, UNMETHSTATE, site.pos_to_end[i], site.prime[i])];
      // deamin_unmethylated = d[1][site.pos_to_end[i]][site.prime[i]];

      if (site.base_compos[i] == 0) {
        // C->C -> (1-e) * (1-d) + (d * e/3)
        M = (1 - seqerror) * (1 - deamin_methylated) + (deamin_methylated * seqerror / 3.0);
        noM = (1 - seqerror) * (1 - deamin_unmethylated) + (deamin_unmethylated * seqerror / 3.0);
      } else if (site.base_compos[i] == 1) {
        // C->T -> (1-e) * d + ((1-d) * e/3)
        M = (1 - seqerror) * (deamin_methylated) + ((1-deamin_methylated) * seqerror / 3.0);
        noM = (1 - seqerror) * deamin_unmethylated + ((1-deamin_unmethylated) * seqerror / 3.0);
      } else {
        // C->[AG] -> (d)*e/3.0 + (1-d) * e/3.0
        M = (seqerror / 3) * deamin_methylated + ((1 - deamin_methylated) * seqerror / 3.0);
        noM = (seqerror / 3) * deamin_unmethylated + ((1 - deamin_unmethylated) * seqerror / 3.0);
      }
      res = (1-maperror) * (d->settings.M * M + (1-d->settings.M)*noM)  + maperror * BASE_FREQ_FLAT_PRIOR;
      ll += std::log(res);
    }
  }

#ifdef skipinternalpos_speak
  std::cerr << "total: " << total_depth << " internal: " << internal_depth << '\n';
#endif

  for(size_t i=0; i< d->nocpg_data->depth; i++){
#ifdef skipinternalpos
    if(nocpg_data.pos_to_end[i]==(settings.max_pos_to_end+1)){
      continue;
    }
#endif
    res = 0;
    seqerror = d->nocpg_data->seqerrors[i];
    maperror = d->nocpg_data->maperrors[i];
    // deamin_methylated = d[0][nocpg_data.pos_to_end[i]][nocpg_data.prime[i]];
    // deamin_unmethylated = d[1][d->nocpg_data->pos_to_end[i]][d->nocpg_data->prime[i]];
    deamin_unmethylated = x[get_param_idx(d->settings.max_pos_to_end, UNMETHSTATE, d->nocpg_data->pos_to_end[i], d->nocpg_data->prime[i])];

    if (d->nocpg_data->base_compos[i] == 0) {
      // C->C -> (1-e) * (1-d) + (d * e/3)
      noM = (1 - seqerror) * (1 - deamin_unmethylated) + (deamin_unmethylated * seqerror / 3.0);
    } else if (d->nocpg_data->base_compos[i] == 1) {
      // C->T -> (1-e) * d + ((1-d) * e/3)
      noM = (1 - seqerror) * deamin_unmethylated       + ((1-deamin_unmethylated) * seqerror / 3.0);
    } else {
      // C->[AG] -> (d)*e/3.0 + (1-d) * e/3.0
      noM = (seqerror / 3) * deamin_unmethylated        + ((1 - deamin_unmethylated) * seqerror / 3.0);
    }
    res = (1-maperror) * noM  + maperror * BASE_FREQ_FLAT_PRIOR;
    ll += std::log(res);
  }
  // uf we are minimizing, replace by -ll
  return ll;
}


std::vector<double> single_array_parameters(const size_t & max_pos_to_end, const mcmc_meth_mat param){
  std::vector<double> res;
  size_t max_length = PRIMES*(max_pos_to_end+2)*PRIMES ;
  res.resize(max_length);
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
        res[get_param_idx(max_pos_to_end, i, p, pr)] = param[i][p][pr];
      }
    }
  }
  return res;
}


// void run_mcmc(const general_settings & settings, const tally_mat_freq & tmf, std::vector<per_site> data){
void run_mcmc(const general_settings & settings, const tally_mat_freq & tmf, std::vector<per_site> data, per_site_nocpg & nocpg_data ){

#ifdef skipinternalpos
  std::cerr << "\t-> Skipping internal data and extending the last estimated parameter for each prime." << '\n';
#endif
  mcmc_meth_mat param_mat = init_mcmc_meth_mat(settings, tmf);
  std::vector<double> param = single_array_parameters(settings.max_pos_to_end, param_mat);
  blabla void_stuff(settings, &data, &nocpg_data);
  size_t n_params = PRIMES * ( settings.max_pos_to_end+2) * PRIMES ;
  nlopt::opt opt(nlopt::LN_BOBYQA, n_params);
  // std::vector<double> lb(n_params, SMALLTOLERANCE), up(n_params, 1-SMALLTOLERANCE);
  std::vector<double> lb(n_params, 0), up(n_params, 1);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(up);

  opt.set_max_objective(objective_func, &void_stuff);
  // if  minimizing remember to set -ll
  // opt.set_min_objective(objective_func, &void_stuff);

  opt.set_xtol_rel(1e-8);
  // opt.set_ftol_rel(1e-8);
  double maxf;
  nlopt::result result = opt.optimize(param, maxf);
  std::cout << "The result is" << std::endl;
  std::cout << result << std::endl;
  std::cout << "After this many iterations: " << void_stuff.iteration << std::endl;
  std::cout << "maxmimal loglikelihood: " << maxf << std::endl;
  for (size_t i=0; i<METHSTATES;i++){
    for (size_t p=0; p<settings.max_pos_to_end+2;p++){
      for (size_t pr=0; pr<PRIMES;pr++){
        if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
          // if((pr==1 && p==0)){
          // we have the center pos as 5 prime 20
          // so no need to print center pos 3 prime 20.
          // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
          continue;
        }
        std::cout << i << " " << p << " " << pr << " " << param[get_param_idx(settings.max_pos_to_end, i, p, pr)] << '\n';
      }
    }
  }

  exit(EXIT_SUCCESS);

  // std::cerr << "CHEATING A BIT" << '\n';
  // std::string fn ="/Users/krishang/Documents/My_PAPERS/MethylationModel/code_mcmc_mle2/dammet_res_good_estimates/20.mcmc.means";
  // std::string fn = "temp/my_deamrates_mapped_10bp_fixedinternal.20.mcmc.means"
  // mcmc_meth_mat d_old = load_deamin_probabilities_temp(settings, fn);
  // mcmc_meth_mat d_new = load_deamin_probabilities_temp(settings, fn);

  // long double ll_old = 0, ll_new=0, alpha=0;
  // ll_old = mcmc_loglike_incl_nocpg(settings, d_old, data, nocpg_data);
  // size_t curr_chain=0, Nchains = 1e7,chainout=500;
  // std::uniform_real_distribution<double> uni_dist(0.0,1.0);
  // size_t total_accepted_steps=0, accepted_steps=0;


  // std::string filename = settings.outbase + ".mcmc";
  // std::ofstream f(filename.c_str());
  // checkfilehandle(f, filename);


  // while (curr_chain < Nchains){
  //   get_new_random_values(settings, d_old, d_new, dist_normal);
  //   // ll_new = mcmc_loglike(settings, d_new, data);
  //   ll_new = mcmc_loglike_incl_nocpg(settings, d_new, data, nocpg_data);

  //   alpha = std::exp(ll_new - ll_old);
  //   if (alpha>1){
  //     d_old = d_new;
  //     ll_old = ll_new;
  //     accepted_steps++;
  //   } else if (uni_dist(gen)<alpha) {
  //     d_old = d_new;
  //     ll_old = ll_new;
  //     accepted_steps++;
  //   } else {
  //     // d_old stays
  //   }
  //   if(curr_chain && curr_chain % chainout == 0){
  //     total_accepted_steps += accepted_steps;
  //     std::cerr << "\t-> MCMC Chain: " << curr_chain << " TotalAcc: " << total_accepted_steps << " Accrate: " << accepted_steps/(double)chainout << " LL: " << ll_old << " Break the chain (:)) with ctrl+c"<< '\r';
  //     accepted_steps=0;
  //     print_mcmc_meth_mat(settings, f, d_old, curr_chain, ll_old);
  //   }
  //   curr_chain++;

  //   // caught a ctrl+c from user. See catch_signal
  //   if (break_MCMC){
  //     std::cerr << "\t-> Caught ctrl+c and will stop the MCMC chain and proceed dammeth" << '\n';
  //     break;
  //   }
  // }
}


std::vector<size_t> load_pos_per_chrom(std::string & filename, const std::string & selected_chrom){
  // 1-based coordinates
  std::vector<size_t> res;
  res.reserve(100000);
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    std::string row, chrom;
    size_t pos;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> chrom >> pos;
      if(chrom == selected_chrom){
        // -1 as it has to zero-based :)
        res.push_back(pos-1);
      }
      ss.clear();
    }
  }
  f.close();
  return res;
}

void filter_ref_db(const std::string & chrom, std::string & filename, char * ref){
  std::vector<size_t> pos = load_pos_per_chrom(filename, chrom);
  for (const auto & val : pos){
    // std::cout << val << " " << ref[val] << '\n';
    ref[val] = 'N';
  }
}

double get_het_prior(const general_settings & settings, const dinucl_pair_of_pair & dinucl_pair){
  double res;
  if(settings.het_rate>=0){
    int difference = ((dinucl_pair.first.first != dinucl_pair.second.first) +
                      (dinucl_pair.first.second != dinucl_pair.second.second));
    if (difference == 0) {
      res = 1 - settings.het_rate - std::pow(settings.het_rate, 2);
    } else if (difference == 1) {
      res = settings.het_rate;
    } else {
      res = std::pow(settings.het_rate, 2);
    }
  } else {  // if het is set to -1, flat prior i.e. 1 -> should be 1/ALL_DINUCL_PAIR_OF_PAIRS.size() but same same
    res = 1;
  }
  return res;
}

void calculate_dinucl_genotypelikelihood(const general_settings & settings, per_site & site, const tally_mat_freq & tmf){
  double res=0, error1_1 = 0, error1_2 = 0, error2_1 = 0, error2_2 = 0;
  double totalres = 0;
  double het_prior = 0;
  for (const auto dinucl_pair : ALL_DINUCL_PAIR_OF_PAIRS) {
    res = 0 ;
    het_prior = get_het_prior(settings, dinucl_pair);

    // readpos prime strand dinucl1.first dinucl1.second dinucl2.first dinucl2.second
    // res[p][pr][strand][r1][r2][s1][s2] = val;
    for (size_t i=0; i<site.bases.size(); i++){
      error1_1 = tmf[site.pos_to_end[i]][site.prime[i]][site.strand[i]]
                    [dinucl_pair.first.first][dinucl_pair.first.second]
                    [site.bases[i].first][dinucl_pair.first.second];
      error1_2 = tmf[site.pos_to_end[i]][site.prime[i]][site.strand[i]]
                    [dinucl_pair.second.first][dinucl_pair.second.second]
                    [site.bases[i].first][dinucl_pair.second.second];
      error2_1 = tmf[site.pos_to_end[i]][site.prime[i]][site.strand[i]]
                    [dinucl_pair.first.first][dinucl_pair.first.second]
                    [dinucl_pair.first.first][site.bases[i].second];
      error2_2 = tmf[site.pos_to_end[i]][site.prime[i]][site.strand[i]]
                    [dinucl_pair.second.first][dinucl_pair.second.second]
                    [dinucl_pair.second.first][site.bases[i].second];

      if (error1_1 < 1e-9) {
        // fprintf(stderr, "these are all the double zeros\n");
        error1_1 = 1e-9;
      }
      if(error1_2<1e-9){
        // fprintf(stderr, "these are all the double zeros\n");
        error1_2=1e-9;
      }
      if(error2_1<1e-9){
        // fprintf(stderr, "these are all the double zeros\n");
        error2_1=1e-9;
      }
      if(error2_2<1e-9){
        // fprintf(stderr, "these are all the double zeros\n");
        error2_2=1e-9;
      }
      res += std::log(((error1_1 * 0.5) + (error1_2 * 0.5)) *
                      ((error2_1 * 0.5) + (error2_2 * 0.5)));
    }
    site.dinucl_genotype_likelihoods.push_back( std::exp(res)*het_prior);
    totalres+=std::exp(res)*het_prior;
    // site.dinucl_genotype_likelihoods.push_back( std::exp(res));
    // totalres+=std::exp(res);
  }
  // divides by the sum of all dinucleotide genotype likelihoods
  for (auto & val: site.dinucl_genotype_likelihoods){
    val = val/totalres;
   }

  // for (size_t i = 1; i < site.dinucl_genotype_likelihoods.size(); i++) {
  //   site.dinucl_genotype_likelihoods[i] = site.dinucl_genotype_likelihoods[i]/totalres;
  // }

  int idx_max = 0;
  for (size_t i = 1; i < site.dinucl_genotype_likelihoods.size(); i++) {
    if (site.dinucl_genotype_likelihoods[i] >
        site.dinucl_genotype_likelihoods[idx_max]) {
      idx_max = i;
    }
  }
  site.dinucl_max_idx = idx_max;
}

mcmc_meth_mat load_deamin_probabilities(const general_settings & settings){
  mcmc_meth_mat res = get_vectors(settings);
  std::string filename = settings.outbase + ".mcmc.means";
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    std::string row;
    size_t methstate, pos, prime;
    double rate;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> methstate >> pos >> prime >> rate;
      ss.clear();
      res[methstate][pos][prime] = rate;
    }
  }
  f.close();
  return res;

}


void dump_deamin_probabilities(const general_settings & settings, mcmc_meth_mat & mat){
  std::string filename = settings.outbase + ".mcmc.means";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    for (size_t i=0; i<METHSTATES;i++){
      for (size_t p=0; p<settings.max_pos_to_end+2;p++){
        for (size_t pr=0; pr<PRIMES;pr++){
          if((pr==1 && p==settings.max_pos_to_end+1) || (pr==1 && p==0)){
          // if((pr==1 && p==0)){
            // we have the center pos as 5 prime 20
            // so no need to print center pos 3 prime 20.
            // no need to estimate pos 0 3 prime as it cannot deaminate. It will alway be a G in a CpG context
            continue;
          }
          f << i << " " << p << " " << pr << " "  << mat[i][p][pr] << '\n';
          // f << methstate << pos << prime << rate;
        }
      }
    }
  }
  f.close();
}

mcmc_meth_mat load_filter_deamin_probabilities(const general_settings & settings, double & burn_in){
  mcmc_meth_mat final_res = get_vectors(settings);
  std::string filename = settings.outbase + ".mcmc";
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    std::string row;
    size_t chain_number=0, prev_chain_number=0;
    size_t methstate, pos, prime;
    double rate;
    std::stringstream ss;
    mcmc_meth_mat res;
    std::vector<mcmc_meth_mat> all_res;
    std::vector<size_t> all_chains;
    res = get_vectors(settings);
    while(getline(f, row)){
      ss.str(row);
      ss >> chain_number >> methstate >> pos >> prime >> rate;
      ss.clear();

      if(chain_number != prev_chain_number && prev_chain_number > 0){
        all_res.push_back(res);
        all_chains.push_back(prev_chain_number);
        res = get_vectors(settings);
      }

      res[methstate][pos][prime] = rate;
      prev_chain_number = chain_number;
    }
    all_res.push_back(res);
    all_chains.push_back(chain_number);


    size_t min_chain = all_chains[all_chains.size()-1] * burn_in;
    size_t chains_used = 0;
    for(size_t i=0; i<all_chains.size(); i++){
      if(all_chains[i] > min_chain){
        chains_used++;
        merge_vectors(settings, final_res, all_res[i]);
      }
    }
    normalize_vectors(settings, final_res, chains_used);

    // here we have to extend the last pos to the internal part of the read;
  }
  f.close();
  return final_res;
}


void update_mle_data(per_mle_run & mle_data, const pre_calc_per_site & data, size_t & idx){
  mle_data.idx_to_include.push_back(idx);
  mle_data.prob_no_deamin_event_meth *= data.prob_no_deamin_event_meth;
  mle_data.prob_no_deamin_event_unmeth *= data.prob_no_deamin_event_unmeth;
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
  for(size_t i=0; i<data.summary.size(); i++){
    mle_data.summary[i] += data.summary[i];
  }
}


struct mle_ll_run {
  double ll, slope, second_der;
};

void mle_loglike(mle_ll_run & d, const double & f, const std::vector<pre_calc_per_site> & data, const per_mle_run & mle_data){
  double ll=0, slope=0, second_der=0;
  double noM, M;
  double mape;
  for (const auto idx : mle_data.idx_to_include){
    for (size_t i=0; i < data[idx].depth; i++){
      noM =  data[idx].pre_noM[i];
      M = data[idx].pre_M[i];
      mape = data[idx].maperrors[i];

      ll += std::log((1-mape) * ((1-f) * noM + f * M) + mape*BASE_FREQ_FLAT_PRIOR);
      // wolfram: derivative  ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // slope += ((1-mape) * (M - noM)) / ((mape-1) * ((f-1) * noM - f * M) + mape*BASE_FREQ_FLAT_PRIOR);
      // beloved wolfram got it right, just modified a few things for readability
      slope += ((1-mape) * (M - noM)) / ((1-mape) * ((1-f) * noM + f * M) + mape*BASE_FREQ_FLAT_PRIOR);
      // https://www.symbolab.com/
      // wolfram: second derivative ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // second derivative:
      // second_der += (std::pow((1-mape),2) * std::pow((M - noM),2)) / std::pow((-M * (mape-1) * f + noM * (mape-1) * (1-f) + mape*BASE_FREQ_FLAT_PRIOR), 2);
    }
  }
  d.ll=ll;
  d.slope=slope;
  // d.second_der = 1.96 * std::sqrt(second_der);
}


void run_mle(const general_settings & settings,
             const std::vector<pre_calc_per_site> & pre_calc_data) {

  for (size_t curr_idx = 0; curr_idx < pre_calc_data.size(); curr_idx++) {

    per_mle_run mle_data(pre_calc_data[curr_idx], curr_idx);
    size_t pos_to_left = (curr_idx > 0) ? curr_idx - 1 : 0;
    size_t pos_to_right = curr_idx + 1;
    bool can_go_left = false, can_go_right = false;
    // keep adding site to mle_data until we fulfill the requirements
    // while ((1 - mle_data.prob_no_deamin_event_meth) < settings.min_prob_deamin_event) {
    // MAX CPG's
    while (mle_data.n_cpgs < 20) {
      can_go_right = false;
      can_go_left = false;

      if (pos_to_right < pre_calc_data.size() - 1) {
        can_go_right = true;
      }

      if (pos_to_left >= 1) {
        can_go_left = true;
      }

      if (!pos_to_right && !pos_to_left) {
        // break ;; cannot extend anymore
        break;

      } else if (pos_to_right && !pos_to_left) {
        // take the one to the right
        update_mle_data(mle_data, pre_calc_data[pos_to_right], pos_to_right);
        pos_to_right++;
      } else if (!pos_to_right && pos_to_left) {
        // take the one to the left
        update_mle_data(mle_data, pre_calc_data[pos_to_left], pos_to_left);
        pos_to_left--;
      } else if (pos_to_right && pos_to_left) {
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
    double ALPHA = 1e-8;
    double tole = 1e-10;
    double BETA = 0.5;

    // here starts the mle estimator.
    double f=0.5, fnew=0.5; // initial guess
    // fprintf(stderr, "startinglikeli: %f\n", oldll);
    size_t niter = 20000000;
    size_t c_iter=0;
    mle_ll_run ll_run, ll_runold;
    mle_loglike(ll_runold, f, pre_calc_data, mle_data);
    fnew = f + ALPHA * ll_runold.slope;
    fprintf(stderr, "startlik: %f f: %f fnew:%f\n", ll_runold.ll, f, fnew);
    double z = ll_runold.slope;
    for (f = 0; f<=1; f+=0.01){
    // for (c_iter = 0 ; c_iter < niter; c_iter++) {
    //   mle_loglike(ll_run, fnew, pre_calc_data, mle_data);
    //   z = BETA * z + ll_run.slope;
    //   fnew = f + ALPHA * z;
      BETA=0;
      mle_loglike(ll_run, f, pre_calc_data, mle_data);
      z = BETA * z + ll_run.slope;
      fnew = f + ALPHA * z;
      fprintf(stdout, "lik=%f diff=%e f=%f fnew=%f slope:%f second_der:%f\n", ll_run.ll, fabs(ll_run.ll - ll_runold.ll), f, fnew, ll_run.slope, ll_run.second_der);
      ll_runold = ll_run;

      // if(c_iter%10000==0 && c_iter)
      //   fprintf(stderr, "[%lu] lik=%f diff=%e f=%f slope:%f\n", c_iter, ll_run.ll, fabs(ll_run.ll - ll_runold.ll),
      //           fnew, ll_run.slope);
      // if (fabs(ll_run.ll - ll_runold.ll) < tole) {
      //   // fprintf(stderr,"breaking\n");
      //   break;
      // }
      // if (fnew>1){
      //   fnew = 1-tole;
      // }
      // if(fnew<0){
      //   fnew = tole;
      // }

      // if (ll_run.ll > ll_runold.ll){
      //   ALPHA *= .9;
      // }

      // f = fnew;
      // ll_runold = ll_run;
    }
    mle_ll_run ll_low, ll_high;
    double low = 0.0000001;
    double high = 0.9999999;
    mle_loglike(ll_low, low, pre_calc_data, mle_data);
    mle_loglike(ll_high, high, pre_calc_data, mle_data);
    std::cout << "Starting EM (contig: " << settings.chrom
              << ") with center position: " << mle_data.curr_pos
              << " ::: N_CpGs: " << mle_data.n_cpgs
              << " ::: Depth: " << mle_data.total_depth
              << " (CG:" << mle_data.summary[0]
              << ",TG:" << mle_data.summary[1]
              << ",{A,G}G:" << mle_data.summary[2]
              << ") ::: Prob_of_deam_meth: " << 1 - mle_data.prob_no_deamin_event_meth
              << " ::: Prob_of_deam_unmeth: " << 1 - mle_data.prob_no_deamin_event_unmeth
              << " ::: Min position: " << mle_data.min_pos
              << " ::: Max position: " << mle_data.max_pos
              << " ::: Distance (bp): " << mle_data.max_pos - mle_data.min_pos
              << " ::: ll: " << ll_run.ll
              << " ::: f: " << f
              << " ::: iterations: " << c_iter
              << " ::: ll0: " << ll_low.ll
              << " ::: ll1: " << ll_high.ll
              << '\n';
    std::cerr << "BREAKING" << '\n';
    break;
  }
}

void filter_per_site_data(const general_settings & settings, const tally_mat_freq & tmf, std::vector<per_site> & data){
    std::vector<per_site> data_temp;
    data_temp.reserve(data.size());
    for (const auto &site : data) {
      // removing sites without depth
      if(site.depth){
        data_temp.push_back(site);
      }
    }
    std::cerr << "\t-> Filtered " << data.size()-data_temp.size() << " sites due to missing data." << '\n';

    // calculate all the dinucleotide genotype likelihoods
    std::cerr << "\t-> Calculating Dinucleotide genotype likelihoods." << '\n';
    for (auto &site : data_temp) {
      calculate_dinucl_genotypelikelihood(settings, site, tmf);
    }
    data.clear();
    // releasing the older data array
    if(true){
      std::cerr << "\t-> Dumping covered " << data_temp.size()  << " dinucls to: " << settings.outbase << ".dinucls" << '\n';
      print_dinucl(settings, data_temp, tmf);
    }
    data.reserve(data_temp.size());
    std::string dinucls_out_fn = settings.outbase + ".dinucls.removed";
    std::ofstream dinucls_out_fh(dinucls_out_fn.c_str());
    checkfilehandle(dinucls_out_fh, dinucls_out_fn);
    std::cerr << "We are not removing any temporarily"<< '\n';
    for (const auto &site : data_temp) {

      if(! check_max_dinucl_CpG(site, MIN_DINUCL_GENOTYPE_PROB) && false){
         dinucls_out_fh << site.position << '\n';
      } else {
        data.push_back(site);
      }
      // std::cerr << "for simple sim" << '\n';
      // data.push_back(site);
    }
    std::cerr << "\t-> " << data_temp.size()-data.size() << " dinucls removed from downstream analyses. Genomic Position dumped to: " << settings.outbase << ".dinucls.removed" << '\n';
}


double calc_prob_obs_base(const double & seqerror, const size_t & base_compos, const double & deamin_rate){
  double res;
  // if(base_compos == 0){  // CC :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
  //   res = ((1-deamin_rate) * ((1-seqerror) * 1) + ((seqerror/3) * seqerror))  + (deamin_rate * ((1-seqerror) * 0) + ((seqerror/3) * seqerror));
  // }else if(base_compos == 1){  // CC :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
  //   res = ((1-deamin_rate) * ((1-seqerror) * 0) + ((seqerror/3) * seqerror))  + (deamin_rate * ((1-seqerror) * 1) + ((seqerror/3) * seqerror));
  // } else { // mutation
  //   res = ((1-deamin_rate) * ((1-seqerror) * 0) + ((seqerror/3) * seqerror))  + (deamin_rate * ((1-seqerror) * 0) + ((seqerror/3) * seqerror));
  // }
  if(base_compos == 0){  // CC :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
    res = ((1-deamin_rate) * (1-seqerror)) + (deamin_rate * seqerror/3) + (seqerror * seqerror/3);
  }else if(base_compos == 1){  // CC :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
    res = (deamin_rate * (1-seqerror)) + ((1-deamin_rate) * seqerror/3) + (seqerror * seqerror/3);
  } else { // mutation
    res = ((1-deamin_rate) * seqerror/3) + (deamin_rate * seqerror/3) + (seqerror * seqerror/3);
  }
  return res;
}


int parse_bam(int argc, char * argv[]) {
  general_settings settings = args_parser(argc, argv);

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

  // go to the correct chromosome
  // https://github.com/gatoravi/bam-parser-tutorial/blob/master/parse_bam.cc
  hts_itr_t *iter = sam_itr_querys(idx, header, settings.chrom.c_str());

  // load ref of that chromosome
  faidx_t *fai = ref_init(settings.reference_fn);

  char *ref = fetch_chrom(fai, settings.chrom);

  char *ref_tally = fetch_chrom(fai, settings.chrom);

  size_t seq_len = faidx_seq_len(fai, settings.chrom.c_str());

  // with true, it works like a charm
  if (!settings.exclude_for_tally_fn.empty()) {
    std::cerr << "\t-> Loading file to mask het and hom alt sites for tally only" << '\n';
    filter_ref_db(settings.chrom, settings.exclude_for_tally_fn, ref_tally);
  }

  if(true){
    std::cerr << "\t-> TEMPORARILY:  REMOVING LOW COVERAGE CPGs from REF" << '\n';
    std::string fn ("temp/meth_low_data_to_be_removed.txt");
    filter_ref_db(settings.chrom, fn, ref);
    filter_ref_db(settings.chrom, fn, ref_tally);
  } else {
    std::cerr << "\t-> TEMPORARILY: NOT REMOVING LOW COVERAGE CPGs from REF" << '\n';
  }

  tally_mat_freq tmf;
  tally_mat tm;
  bool count_table_exists=check_file_exists(settings.outbase + ".tallycounts");

  if (count_table_exists) {
    std::cerr << "\t-> Loading Count file: " << settings.outbase << ".tallycounts" << '\n';
    tmf = read_count_file(settings);
  } else {
    std::cerr << "\t-> Init Count matrix" << '\n';
    tm = init_tallymat(settings.max_pos_to_end);
  }

  const keeplist_map cpg_map = get_cpg_genomic_pos(ref, seq_len);
  std::cerr << "\t-> " << cpg_map.size() << " CpG's in chrom: " << settings.chrom << '\n';

  // Initiate the alignment record
  std::vector<per_site> data;
  data.resize(cpg_map.size());

  size_t nocpg_to_include = 2e6;
  std::cerr << "\t-> Loading " << nocpg_to_include << " of C's outside a CpG context" << '\n';
  per_site_nocpg nocpg_data(nocpg_to_include);

  bam1_t *rd = bam_init1();
  int reads = 1e6;
  int counter = 0;
  alignment_data d;
  size_t trashed=0;
  while (sam_itr_next(in, iter, rd) >= 0) {
    if (counter % reads == 0) {
      std::cerr << "\t-> " << counter << " reads processed and " << trashed << " discarded. " << '\r';
    }

    if (rd->core.qual < settings.minmapQ || (rd->core.flag & settings.flags_off)){
      trashed++;
      continue;
    }

    d = align_read(rd, ref);
    if(! count_table_exists ){
      tally_aligned_data(settings, d, tm, ref_tally);
    }
    add_aligned_data(settings, d, data, cpg_map, nocpg_data, ref_tally);
    // add_aligned_data(settings, d, data, cpg_map);
    counter++;
  }

  if (! count_table_exists ) {
    std::cerr << "\t-> Dumping count file: " << settings.outbase << ".tallycounts" << '\n';
    dump_count_file(settings, tm);
    tmf = read_count_file(settings);
  }
  std::cerr << "\t-> Processed: " << counter+trashed  << ". Reads filtered: " << trashed << '\n';
  filter_per_site_data(settings, tmf, data);
  mcmc_meth_mat mcmc_mat;

  // start MCMC
  // if(check_file_exists(settings.outbase+".mcmc.means") && false){
  // std::string fn;
  // mcmc_meth_mat d_old;
  // double ll_old;
  // std::cerr.precision(30);
  // fn = "/Users/krishang/Documents/My_PAPERS/MethylationModel/code_mcmc_mle2/20.mcmc.means.perfect";
  //  d_old = load_deamin_probabilities_temp(settings, fn);
  // ll_old = mcmc_loglike_incl_nocpg(settings, d_old, data, nocpg_data);
  // std::cerr << "PERFECT: " << ll_old << '\n';
  // fn ="/Users/krishang/Documents/My_PAPERS/MethylationModel/code_mcmc_mle2/dammet_res_good_estimates_fiddled/20.mcmc.means";
  // d_old = load_deamin_probabilities_temp(settings, fn);
  // ll_old = mcmc_loglike_incl_nocpg(settings, d_old, data, nocpg_data);
  // std::cerr << "FIDDLED: " << ll_old << '\n';
  // fn ="/Users/krishang/Documents/My_PAPERS/MethylationModel/code_mcmc_mle2/dammet_res_good_estimates/20.mcmc.means";
  // d_old = load_deamin_probabilities_temp(settings, fn);
  // ll_old = mcmc_loglike_incl_nocpg(settings, d_old, data, nocpg_data);
  // std::cerr << "ESTIMATED: " << ll_old << '\n';

  // exit(EXIT_SUCCESS);


  if(check_file_exists(settings.outbase+".mcmc.means")){
    std::cerr << "\t-> loading parameters from " << settings.outbase << ".mcmc.means" << '\n';
    mcmc_mat = load_deamin_probabilities(settings);
    // std::string fn="20.mcmc.means.perfect.fixed_internal";
    // std::string fn="20.mcmc.means.perfect.fixed_internal_plus010";
    // std::string fn="/Users/krishang/Documents/My_PAPERS/MethylationModel/code_mcmc_mle2/deam_rates/perfect.txt";
    std::string fn="/Users/krishang/Documents/My_PAPERS/MethylationModel/code_mcmc_mle2/deam_rates_20/perfect.txt";
    std::cerr << "\t-> loading PERFECT parameters from " << fn << '\n';
    mcmc_mat = load_deamin_probabilities_temp(settings, fn);
    mcmc_mat[0][settings.max_pos_to_end+1][0] = mcmc_mat[0][settings.max_pos_to_end][0];
    mcmc_mat[1][settings.max_pos_to_end+1][0] = mcmc_mat[1][settings.max_pos_to_end][0];

  } else if(check_file_exists(settings.outbase+".mcmc")){
    // removing 25% burnin and using remaining for the deamin_probabilities
    std::cerr << "\t-> Filtering for burn-in and calculating mean deamination probability." << '\n';
    mcmc_mat = load_filter_deamin_probabilities(settings, BURN_IN_FRAC);
    dump_deamin_probabilities(settings, mcmc_mat);
  } else {
    std::cerr << "\t-> Starting MCMC." << '\n';
    run_mcmc(settings, tmf, data, nocpg_data);
    // run_mcmc(settings, tmf, data);
    std::cerr << "\t-> Filtering for burn-in and calculating mean deamination probability." << '\n';
    // removing 25% burnin and using remaining for the deamin_probabilities
    mcmc_mat = load_filter_deamin_probabilities(settings, BURN_IN_FRAC);
    dump_deamin_probabilities(settings, mcmc_mat);
  }
  std::vector<pre_calc_per_site> mle_data;
  // std::cerr << "\t-> does not run MLE." << '\n';
  if (true) {
    double deamin_methylated, deamin_unmethylated;
    double noM, M;
    std::cerr << "\t-> Precalculating many things for ML." << '\n';
    mle_data.reserve(data.size());
    for (const auto &site : data) {
      pre_calc_per_site d(site);
      for (size_t i = 0; i < site.depth; i++) {
        deamin_methylated = mcmc_mat[0][site.pos_to_end[i]][site.prime[i]];
        deamin_unmethylated = mcmc_mat[1][site.pos_to_end[i]][site.prime[i]];
        noM = calc_prob_obs_base(site.seqerrors[i], site.base_compos[i], deamin_unmethylated);
        M = calc_prob_obs_base(site.seqerrors[i], site.base_compos[i], deamin_methylated);

        d.pre_M.push_back(M);
        d.pre_noM.push_back(noM);

        d.maperrors.push_back(site.maperrors[i]);

        d.prob_no_deamin_event_meth *= 1-deamin_methylated;
        d.prob_no_deamin_event_unmeth *= 1-deamin_unmethylated;
        d.summary[site.base_compos[i]]++;

        if(site.position<=62791){
          std::cout << site.base_compos[i] << " " << site.seqerrors[i] << " pos&prime " << site.pos_to_end[i] << " " << site.prime[i] << " deam&MnoM: " <<  deamin_methylated << " " << deamin_unmethylated << " M&noM: " << M << " " << noM << '\n';
        }
        // if(site.base_compos[i]!=1){
        //   d.pre_M.push_back(M);
        //   d.pre_noM.push_back(noM);

        //   d.maperrors.push_back(site.maperrors[i]);

        //   d.prob_no_deamin_event_meth *= 1-deamin_methylated;
        //   d.prob_no_deamin_event_unmeth *= 1-deamin_unmethylated;
        //   d.summary[site.base_compos[i]]++;

        //   if(site.position<=62791){
        //     std::cout << site.base_compos[i] << " " << site.seqerrors[i] << " pos&prime " << site.pos_to_end[i] << " " << site.prime[i] << " deam&MnoM: " <<  deamin_methylated << " " << deamin_unmethylated << " M&noM: " << M << " " << noM << '\n';
        //   }

        // } else {
        //   d.pre_M.push_back(1);
        //   d.pre_noM.push_back(1);
        //   d.maperrors.push_back(1);

        // }

      }

      mle_data.push_back(d);
    }

    data.clear();
    std::cerr << "\t-> Starting MLE." << '\n';
    run_mle(settings, mle_data);
  }

  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(rd);
  bam_hdr_destroy(header);
  sam_close(in);

  return 0;
}

int main(int argc, char *argv[]) {
  return parse_bam(argc, argv);

}
