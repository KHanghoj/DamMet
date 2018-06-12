#include "deammeth.h"
#include "nlopt.hpp"
#include <set>
#include <ctime>                                                                                                                                                                                             

std::string UNKNOWN_RG("UNKNOWN");
std::string ALL_RG("ALL_DEAMMETH");
size_t ALL_DEAMMETH_RG_IDX = 0;
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

double phred_to_double(int & phred){
  return std::pow(10, -((double)phred/10.0));
}

std::vector<double> phred_to_double_converter(){
  std::vector<double> res;
  for (int phred=0; phred<=255; phred++){
    res.push_back(phred_to_double(phred));
  }
  return res;
}

std::vector<double> PHRED_TO_PROB_CONVERTER = phred_to_double_converter();

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


// tmf: // readpos, prime, strand, refdinucl1, refdinucl2, sampledinucl1,sampledinucl2

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

const std::vector<int> get_cpg_bool(const char * ref, const size_t & seq_len){
  std::vector<int> res(seq_len, 0);
  for (size_t i=0; i<seq_len-1; i++){
    if (refToInt[(int)ref[i]] == 1 && refToInt[(int)ref[i+1]] == 2){
      res[i] = 1;
      res[i+1] = 1;
    }
  }
  return res;
}


void dump_count_file(general_settings & settings, tally_mat & tm, std::string & rg){
  std::string filename = settings.outbase + "." + rg + ".tallycounts";
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

tally_mat_freq read_count_file(general_settings & settings, std::string & rg){
  tally_mat temp = init_tallymat(settings.max_pos_to_end);
  tally_mat_freq res = init_tallymat_freq(settings.max_pos_to_end);
  std::string filename = settings.outbase + "." + rg + ".tallycounts";
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
            size_t totalsum = 0;
            for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
              for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                if (temp[p][pr][strand][r1][r2][s1][s2]) {
                  totalsum += temp[p][pr][strand][r1][r2][s1][s2];
                }
              }
            }
            for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
              for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                if (temp[p][pr][strand][r1][r2][s1][s2] && totalsum) {
                  res[p][pr][strand][r1][r2][s1][s2] = (temp[p][pr][strand][r1][r2][s1][s2] / (double)totalsum);
                }

                if (res[p][pr][strand][r1][r2][s1][s2] < SMALLTOLERANCE){
                  res[p][pr][strand][r1][r2][s1][s2] = SMALLTOLERANCE;
                }
                if(res[p][pr][strand][r1][r2][s1][s2] > 1-SMALLTOLERANCE){
                  res[p][pr][strand][r1][r2][s1][s2] = 1-SMALLTOLERANCE;
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

// Printing dinucl_data
void print_dinucl(const general_settings &settings,
                  const std::vector<per_site> &data) {

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
}

bool nucleotide_missing(const int & r_base1, const int & r_base2, const int & s_base1, const int & s_base2){
  return (r_base1==4 || r_base2==4 || s_base1==4 || s_base2==4 || r_base1==9 || r_base2==9 || s_base1==9 || s_base2==9);
}

bool base_is_indel(const int & base1, const int & base2){
  return (base1==9 || base2==9);
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
          // d.t_ref.push_back(4);
	  // 9 reflects an indel
          d.t_ref.push_back(9);

          d.t_positions.push_back(0);

          seq_pos++; // <- important, must be after macro
        }
        wpos++;
      } else { // this is the deletion part
        for(int ii=0;ii<opLen;ii++){
          // d.t_seq.push_back(4);
	  // 9 reflects an indel
          d.t_seq.push_back(9);
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
}

bool cpg(const int & base1, const int & base2){
  return(base1==1 && base2==2);
}

void add_aligned_data(const general_settings &settings,
                      const alignment_data &d,
                      const keeplist_map &cpg_map,
                      const char *ref,
                      std::vector<per_site> &data,
                      tally_mat &tm,
                      per_site_nocpg &nocpg_data,
		      my_cov_rg & cov_rg,
		      size_t & rg_idx,
		      size_t & ncycles) {
  int b1, b2, r_base1, r_base2, s_base1, s_base2;
  size_t dist_5p, dist_3p;
  size_t pos, prime;

  per_site *p_per_site;
  size_t data_idx=0;
  double seqerror = 0;
  double maperror = PHRED_TO_PROB_CONVERTER[d.mapQ];
  int base_composition; // 0,1,2
  bool include_in_cpg = false;
  bool include_in_tally = false;
  for (size_t i = 0; i < d.t_seq.size()-1; i++){
    include_in_cpg = true;
    include_in_tally = true;

    if(d.strand){ // negative strand
      b1 = d.t_seq.size()-i-2;
      b2 = d.t_seq.size()-i-1; // closest to 5 prime
      get_bases(d, b1, b2, r_base1, r_base2, s_base1, s_base2);

      if(nucleotide_missing(r_base1, r_base2, s_base1, s_base2) ||  d.t_qs[b1] < settings.minbaseQ || d.t_qs[b2] < settings.minbaseQ ){
        continue;
      }

      // check if left side has an indel
      if(b1-1 >=0 && base_is_indel(d.t_ref[b1-1], d.t_seq[b1-1])){
      	continue;
      }
      
      // check if right side has an indel
      if(b2+1 < d.t_seq.size()-1 && base_is_indel(d.t_ref[b2+1], d.t_seq[b2+1])){
      	continue;	
      }

      dist_5p = d.t_isop[b2];
      dist_3p = d.t_posi[b2];
      seqerror = PHRED_TO_PROB_CONVERTER[d.t_qs[b2]];
      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

      // discard 3prime data, if read == number of cycles.
      if (prime==1 && d.n_nucleotides==ncycles){
	// instead of discarding, we should just use the postions relative to the 5 prime. 
	// if 40 bp from 5prime and 30bp from 3 prime, we will just use the 40bp 5 prime.
	get_pos_and_prime(settings, dist_5p*10, dist_5p, pos, prime);
	// continue;
      }

      if(refToInt[(int)ref[d.t_positions[b1]]]==4 || refToInt[(int)ref[d.t_positions[b2]]]==4){
        include_in_tally=false;
      }


      if (s_base2 == 2){
        base_composition = 0;  // no error
      } else if (s_base2 == 0){
        base_composition = 1;  // deamin or error
      } else {
        base_composition = 2; // error (could also be deamin then error)
      }

      if(no_ref_cpg(r_base1, r_base2)){
        if(refToInt[(int)ref[d.t_positions[b2]]]==2 && cov_rg.nocpg[rg_idx] < cov_rg.cpg[rg_idx]){
          nocpg_data.pos_to_end.push_back(pos);
          nocpg_data.prime.push_back(prime);
          nocpg_data.base_compos.push_back(base_composition);
          nocpg_data.seqerrors.push_back(seqerror);
          nocpg_data.maperrors.push_back(maperror);
	  nocpg_data.rgs.push_back(rg_idx);
	  nocpg_data.depth++;
	  cov_rg.nocpg[rg_idx]++;
        }
        include_in_cpg = false;
      }
    } else { // positive strand
      b1 = i; // closest to 5 prime
      b2 = i+1;
      get_bases(d, b1, b2, r_base1, r_base2, s_base1, s_base2);

      if(nucleotide_missing(r_base1, r_base2, s_base1, s_base2) ||  d.t_qs[b1] < settings.minbaseQ || d.t_qs[b2] < settings.minbaseQ ){
        continue;
      }

      // check if left side has an indel
      if(b1-1 >=0 && base_is_indel(d.t_ref[b1-1], d.t_seq[b1-1])){
      	continue;
      }
      
      // check if right side has an indel
      if(b2+1 < d.t_seq.size()-1 && base_is_indel(d.t_ref[b2+1], d.t_seq[b2+1])){
      	continue;	
      }

      dist_5p = d.t_posi[b1];
      dist_3p = d.t_isop[b1];
      seqerror = PHRED_TO_PROB_CONVERTER[d.t_qs[b1]];
      get_pos_and_prime(settings, dist_3p, dist_5p, pos, prime);

      // discard 3prime data, if read == number of cycles.
      if (prime==1 && d.n_nucleotides==ncycles){
	// instead of discarding, we should just use the postions relative to the 5 prime. 
	// if 40 bp from 5prime and 30bp from 3 prime, we will just use the 40bp 5 prime.
	get_pos_and_prime(settings, dist_5p*10, dist_5p, pos, prime);
	// continue;
      }

      if(refToInt[(int)ref[d.t_positions[b1]]]==4 || refToInt[(int)ref[d.t_positions[b2]]]==4){
        include_in_tally=false;
      }

      if (s_base1 == 1){
        base_composition = 0;  // no error
      } else if (s_base1 == 3){
        base_composition = 1;  // deamin or error
      } else {
        base_composition = 2; // error (could also be deamin then error)
      }

      if(no_ref_cpg(r_base1, r_base2)){
        if(refToInt[(int)ref[d.t_positions[b1]]]==1 && cov_rg.nocpg[rg_idx] < cov_rg.cpg[rg_idx]){
          nocpg_data.pos_to_end.push_back(pos);
          nocpg_data.prime.push_back(prime);
          nocpg_data.base_compos.push_back(base_composition);
          nocpg_data.seqerrors.push_back(seqerror);
          nocpg_data.maperrors.push_back(maperror);
	  nocpg_data.rgs.push_back(rg_idx);
	  nocpg_data.depth++;
	  cov_rg.nocpg[rg_idx]++;
        }
        include_in_cpg = false;
      }

    }
    if(include_in_tally){
      tm[pos][prime][d.strand][r_base1][r_base2][s_base1][s_base2]++;
    }
    if(include_in_cpg){
      data_idx = cpg_map.at(d.t_positions[b1]);
      p_per_site = &data[data_idx];
      p_per_site->position  = d.t_positions[b1];

      p_per_site->depth++;
      p_per_site->strand.push_back(d.strand);

      p_per_site->pos_to_end.push_back(pos);
      p_per_site->prime.push_back(prime);

      p_per_site->base_compos.push_back(base_composition);
      p_per_site->seqerrors.push_back(seqerror);
      p_per_site->maperrors.push_back(maperror);
      p_per_site->rgs.push_back(rg_idx);
      cov_rg.cpg[rg_idx]++;

      p_per_site->refs.push_back(std::make_pair(r_base1, r_base2));
      p_per_site->bases.push_back(std::make_pair(s_base1, s_base2));
      p_per_site->quals.push_back(std::make_pair(d.t_qs[b1], d.t_qs[b2]));
    }
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

struct F_void {
  general_settings settings;
  std::vector<pre_calc_per_site> * data;
  per_mle_run * mle_data;
  size_t iteration;
  F_void(const general_settings & s, std::vector<pre_calc_per_site> * d, per_mle_run * m){
    settings = s;
    data = d;
    mle_data = m;
    iteration = 0;
  }
};


struct deamrates_void {
  general_settings settings;
  std::vector<per_site> * data;
  per_site_nocpg * nocpg_data;
  size_t iteration, rg_idx;
  std::vector< std::vector< size_t > > * cpg_idx;
  std::vector< size_t > * nocpg_idx;
  deamrates_void(const general_settings & s, std::vector<per_site> * d, per_site_nocpg * nocpg_d, size_t & r_idx, std::vector< std::vector< size_t > > * c_idx, std::vector< size_t > * nc_idx){
    settings = s;
    data = d;
    nocpg_data = nocpg_d;
    rg_idx = r_idx;
    cpg_idx = c_idx;
    nocpg_idx = nc_idx;
    iteration = 0;
  }
};

// struct deamrates_void {
//   general_settings settings;
//   std::vector<per_site> * data;
//   per_site_nocpg * nocpg_data;
//   size_t iteration, rg_idx;
//   deamrates_void(const general_settings & s, std::vector<per_site> * d, per_site_nocpg * nocpg_d, size_t & idx){
//     settings = s;
//     data = d;
//     nocpg_data = nocpg_d;
//     rg_idx = idx;
//     iteration = 0;
//   }
// };


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
  per_site * p_per_site;
  for (size_t site=0; site<d->cpg_idx->size(); site++){
    // http://tytel.org/blog/2015/01/19/unexpected-optimization-00/
    // p_per_site = &((*d->data)[site]);
    p_per_site = &d->data->at(site);
    // if((*d->cpg_idx)[site].size()==0){
    //   continue;
    // }
    for(const auto & i : d->cpg_idx->at(site)){
    // for(const auto & i : (*d->cpg_idx)[site]){

      idx_param_meth = get_param_idx(d->settings.max_pos_to_end, METHSTATE, p_per_site->pos_to_end[i], p_per_site->prime[i]);
      idx_param_unmeth = get_param_idx(d->settings.max_pos_to_end, UNMETHSTATE, p_per_site->pos_to_end[i], p_per_site->prime[i]);
      res=0;
      seqerror = p_per_site->seqerrors[i];
      maperror = p_per_site->maperrors[i];
      deamin_methylated = x[idx_param_meth];
      deamin_unmethylated = x[idx_param_unmeth];
      map_and_prior =  maperror * BASE_FREQ_FLAT_PRIOR;

      if (p_per_site->base_compos[i] == 0) {
        // C->C -> (1-e) * (1-d) + (d * e/3)
        M = (1 - seqerror) * (1 - deamin_methylated) + (deamin_methylated * seqerror / 3.0);
        noM = (1 - seqerror) * (1 - deamin_unmethylated) + (deamin_unmethylated * seqerror / 3.0);
        res = (1-maperror) * (d->settings.M * M + (1-d->settings.M)*noM) + map_and_prior;

        if (!grad.empty()) {
          grad[idx_param_meth] += (d->settings.M*(1.0-maperror)*( 4.0*seqerror / 3.0 - 1) + map_and_prior) / res;
          grad[idx_param_unmeth] += ((1-d->settings.M) * (1.0-maperror)*( 4.0*seqerror / 3.0 - 1) + map_and_prior) /res;
        }

      } else if (p_per_site->base_compos[i] == 1) {
        // C->T -> (1-e) * d + ((1-d) * e/3)
        M = (1 - seqerror) * (deamin_methylated) + ((1-deamin_methylated) * seqerror / 3.0);
        noM = (1 - seqerror) * deamin_unmethylated + ((1-deamin_unmethylated) * seqerror / 3.0);
        res = (1-maperror) * (d->settings.M * M + (1-d->settings.M)*noM)  + map_and_prior;
        if (!grad.empty()) {
          grad[idx_param_meth] += (d->settings.M * (1.0-maperror)*( -4.0*seqerror / 3.0 + 1) + map_and_prior) / res;
          grad[idx_param_unmeth] += ((1-d->settings.M) * (1.0-maperror)*( -4.0*seqerror / 3.0 + 1) + map_and_prior) / res;
        }
      } else {
        // C->[AG] -> (d)*e/3.0 + (1-d) * e/3.0
        M = (seqerror / 3) * deamin_methylated + ((1 - deamin_methylated) * seqerror / 3.0);
        noM = (seqerror / 3) * deamin_unmethylated + ((1 - deamin_unmethylated) * seqerror / 3.0);
        res = (1-maperror) * (d->settings.M * M + (1-d->settings.M)*noM)  + map_and_prior;

        // diff is zero
        // grad[idx_param_unmeth] += 0;
        // grad[idx_param_meth] += 0;

      }
      ll += std::log(res);

    }
  }

  // if (!grad.empty()) {
  //   std::cerr << "before" << '\n';
  //   print_single_array_parameters(d->settings.max_pos_to_end, grad, std::cerr);
  //   std::cerr << "after" << '\n';
  // }
  for(const auto & i : (*d->nocpg_idx)){

    idx_param_unmeth = get_param_idx(d->settings.max_pos_to_end, UNMETHSTATE, d->nocpg_data->pos_to_end[i], d->nocpg_data->prime[i]);
    // std::cerr << "got tis far " << d->nocpg_data->depth << " " << d->settings.max_pos_to_end << " " << i << " " << d->nocpg_data->pos_to_end[i]  << " " << d->nocpg_data->prime[i] << " " << idx_param_unmeth<< " " << x.size() <<  " " << d->nocpg_data->seqerrors.size() << " " << d->nocpg_data->maperrors.size() << '\n';
    res = 0;
    seqerror = d->nocpg_data->seqerrors[i];
    maperror = d->nocpg_data->maperrors[i];
    deamin_unmethylated = x[idx_param_unmeth];
    map_and_prior =  maperror * BASE_FREQ_FLAT_PRIOR;
    
    if (d->nocpg_data->base_compos[i] == 0) {
      // C->C -> (1-e) * (1-d) + (d * e/3)
      noM = (1 - seqerror) * (1 - deamin_unmethylated) + (deamin_unmethylated * seqerror / 3.0);
      res = (1-maperror) * noM  + map_and_prior;

      if (!grad.empty()) {
        grad[idx_param_unmeth] += ((1-maperror) * ( 4.0*seqerror / 3.0 - 1) + map_and_prior) / res;
      }

    } else if (d->nocpg_data->base_compos[i] == 1) {
      // C->T -> (1-e) * d + ((1-d) * e/3)
      noM = (1 - seqerror) * deamin_unmethylated + ((1-deamin_unmethylated) * seqerror / 3.0);
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
  //   print_single_array_parameters(d->settings.max_pos_to_end, grad, std::cerr);
  // }
  // return -ll;
  return ll;
}


std::vector<double> single_array_parameters(const size_t & max_pos_to_end, const tally_mat_freq & tmf){
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
          res[get_param_idx(max_pos_to_end, i, p, pr)] = tmf[p][pr][0][1][2][3][2];
        } else {
          res[get_param_idx(max_pos_to_end, i, p, pr)] = tmf[p][pr][0][1][1][3][1];
        }
      }
    }
  }
  return res;
}

// -109750 with BOBYQA 10 minutes
// -109747 with MMA (requires derivatives) 20 seconds wuhu
// void run_deamrates_optim(const general_settings & settings, const tally_mat_freq & tmf, std::vector<per_site> data){
void run_deamrates_optim(const general_settings & settings, const tally_mat_freq & tmf, std::vector<per_site> & data, per_site_nocpg & nocpg_data , std::string & rg, size_t & rg_idx, my_cov_rg & cov_rg){
  time_t start_time_optim,end_time_optim;
  time(&start_time_optim);
  std::string filename = settings.outbase + "." + rg + ".deamrates";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  std::vector<double> param = single_array_parameters(settings.max_pos_to_end, tmf);

  std::vector<size_t> nocpg_idx;
  nocpg_idx.reserve(cov_rg.nocpg[rg_idx]);
  for (size_t i=0; i<nocpg_data.depth; i++){
    if(nocpg_data.rgs[i] == rg_idx){
      nocpg_idx.push_back(i);
    }
  }

  std::vector< std::vector< size_t > > cpg_idx;
  cpg_idx.resize(data.size());
  for (size_t site=0; site<data.size(); site++){
    for(size_t i=0; i<data[site].depth; i++){
      if(data[site].rgs[i]==rg_idx){
	cpg_idx[site].push_back(i);
      }
    }
  }

  
  // deamrates_void void_stuff(settings, &data, &nocpg_data, rg_idx);
  deamrates_void void_stuff(settings, &data, &nocpg_data, rg_idx, &cpg_idx, &nocpg_idx);
  // print_single_array_parameters(settings.max_pos_to_end, param, std::cerr);
  std::vector<double> t; // just a dummy
  f << "##Initial (param are freq from tallycounts file) llh: " << objective_func_deamrates(param, t, &void_stuff) << std::endl;
  size_t n_params = PRIMES * ( settings.max_pos_to_end+2) * PRIMES ;
  // nlopt::opt opt(nlopt::LN_BOBYQA, n_params);
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
  std::cerr << "\t-> Iteration: " << void_stuff.iteration << " llh: " << minf  << " in " << difftime(end_time_optim, start_time_optim) << " seconds" << '\n';
  f << "##Optim return code: " << result << std::endl;
  f << "##Iterations: " << void_stuff.iteration << std::endl;
  f << "##llh: " << minf << std::endl;

  print_single_array_parameters(settings.max_pos_to_end, param, f);
  //   print_single_array_parameters(settings.max_pos_to_end, param, std::cout);
  f.close();
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
    } else { // double het
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

std::vector<double> load_deamrates(const general_settings & settings, std::string & rg){
  size_t max_length = PRIMES*(settings.max_pos_to_end+2)*PRIMES ;
  std::vector<double> res(max_length, SMALLTOLERANCE);
  std::string filename = settings.outbase + "." + rg + ".deamrates";
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
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

std::vector<double> load_deamrates_f(const general_settings & settings){
  size_t max_length = PRIMES*(settings.max_pos_to_end+2)*PRIMES ;
  std::vector<double> res(max_length, SMALLTOLERANCE);
  std::string filename = settings.deamrates_filename;
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
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

double objective_func_F(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
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

// this is to crazy higher. 95% confidence intervals of 7.6658. no good at all
double objective_func_F_second_deriv(const double & f, std::vector<pre_calc_per_site> & data, per_mle_run &mle_data){
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


void run_mle(const general_settings & settings,
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

  // nlopt::opt opt(nlopt::LN_BOBYQA, n_params);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  std::vector<double> lb(1, SMALLTOLERANCE), up(1, 1-SMALLTOLERANCE);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(up);
  opt.set_xtol_abs(0);
  opt.set_ftol_abs(1e-30);
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
    F_void void_stuff(settings, &pre_calc_data, &mle_data);
    std::vector<double> param (1, 0.5);
    nlopt::result result;
    double second_der, error;
    size_t iterations;
      
    if(! same){
      opt.set_max_objective(objective_func_F, &void_stuff);
      result = opt.optimize(param, minf);
      second_der = objective_func_F_second_deriv(param[0], pre_calc_data, mle_data);
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

    f << "contig: " << settings.chrom
      << " Center_pos: " << mle_data.curr_pos
      << " N_CpGs: " << mle_data.n_cpgs
      << " Depth: " << mle_data.total_depth
      << " Stats: " << mle_data.summary[0] << "," << mle_data.summary[1] << "," << mle_data.summary[2]
      << " Prob(deam_m): " << 1 - mle_data.prob_no_deamin_event_meth
      << " Incl_pos: ";
    // std::sort(mle_data.positions.begin(), mle_data.positions.end());
    f << mle_data.positions[0];
    for (auto i=mle_data.positions.begin()+1; i!=mle_data.positions.end(); i++){
      f << ","<< *i;
    }
    f << " Distance: " << mle_data.max_pos - mle_data.min_pos
      << " ll: " << minf
      << " f: " << param[0]
      << " f(95%conf): " << param[0]-error  << "," <<  param[0]+error
      << " iterations: " << iterations
      << " return_code: " << result
      << '\n';
  }
  f << std::flush;
  f.close();
}

void filter_per_site_data(const general_settings & settings, const tally_mat_freq & tmf, std::vector<per_site> & data, std::string & rg){
    std::vector<per_site> data_temp;
    data_temp.reserve(data.size());
    for (const auto &site : data) {
      // removing sites without depth
      if(site.depth){
        data_temp.push_back(site);
      }
    }
    
    std::cerr << "\t-> Filtered RG: " << rg << " " << data.size()-data_temp.size() << " sites due to missing data." << '\n';
    // calculate all the dinucleotide genotype likelihoods
    std::cerr << "\t-> TEMPORARILY:  NOT CALCULATING DINUCL GENOTYPE LIKELIHHODS. Only removing non-covered cpgs" << '\n';
    data.clear();
    data.reserve(data_temp.size());
    for (const auto &site : data_temp) {
      data.push_back(site);
    }

    
    if(false){
      std::cerr << "\t-> Calculating Dinucleotide genotype likelihoods." << '\n';
      for (auto &site : data_temp) {
	calculate_dinucl_genotypelikelihood(settings, site, tmf);
      }
      data.clear();
      // releasing the older data array
      if(true){
	std::cerr << "\t-> Dumping covered " << data_temp.size()  << " dinucls to: " << settings.outbase << ".dinucls" << '\n';
	print_dinucl(settings, data_temp);
      }
      data.reserve(data_temp.size());
      std::string dinucls_out_fn = settings.outbase + ".dinucls.removed";
      std::ofstream dinucls_out_fh(dinucls_out_fn.c_str());
      checkfilehandle(dinucls_out_fh, dinucls_out_fn);
      std::cerr << "TEMPORARILY: We are not removing any temporarily"<< '\n';
      // std::cerr << "YES WE ARE TEMPORARILY"<< '\n';
      for (const auto &site : data_temp) {
	
	if(! check_max_dinucl_CpG(site, MIN_DINUCL_GENOTYPE_PROB) && false){
	  // if(! check_max_dinucl_CpG(site, MIN_DINUCL_GENOTYPE_PROB)){
	  dinucls_out_fh << site.position << '\n';
	} else {
	  data.push_back(site);
	}
      }
      std::cerr << "\t-> " << data_temp.size()-data.size() << " dinucls removed from downstream analyses. Genomic Position dumped to: " << settings.outbase << ".dinucls.removed" << '\n';
    }
    
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

void dump_options(const general_settings & settings){
  std::string filename = settings.outbase+".args";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  f<<settings.all_options;
  f.close();
}

int parse_bam(int argc, char * argv[]) {
  time_t start_time, end_time;
  time(&start_time);
  general_settings settings = args_parser(argc, argv);
  dump_options(settings);
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
  if(!settings.deamrates_f_rg.empty() && check_file_exists(settings.deamrates_f_rg)){
    rg_split = true;
    std::ifstream f (settings.deamrates_f_rg.c_str());
    checkfilehandle(f, settings.deamrates_f_rg);
    std::string row, rg, cycle_string;
    size_t cycle;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> rg >> cycle_string;
      if(cycle_string.empty()){
	cycle=std::numeric_limits<size_t>::max();
	std::cerr << "\t-> RG: " << rg  << ". setting sequencing cycles: " << cycle << '\n';
      } else {
	cycle=std::stoi(cycle_string);
	std::cerr << "\t-> RG: " << rg  << ". sequencing cycles: " << cycle << '\n';
      }
      rgs.push_back(rg);
      cycles.push_back(cycle);
      ss.clear();
    }
    f.close();
  } else {
    rg_split = false;
    rgs.push_back(ALL_RG);
    if(settings.cycles!=std::numeric_limits<size_t>::max()){
      cycles.push_back(settings.cycles);
    } else {
      cycles.push_back(std::numeric_limits<size_t>::max());
    }
    std::cerr << "\t-> Merging all reads into a single read group named: " << rgs[0] << " and sequencing cycles: " << cycles[0]  << '\n';
  }

  // go to the correct chromosome
  // https://github.com/gatoravi/bam-parser-tutorial/blob/master/parse_bam.cc
  hts_itr_t *iter = sam_itr_querys(idx, header, settings.chrom.c_str());

  // load ref of that chromosome
  faidx_t *fai = ref_init(settings.reference_fn);

  char *ref = fetch_chrom(fai, settings.chrom);

  size_t seq_len = faidx_seq_len(fai, settings.chrom.c_str());

  // with true, it works like a charm
  if (!settings.exclude_sites_fn.empty()) {
    std::cerr << "\t-> Loading file to mask genomic sites sites" << '\n';
    filter_ref_db(settings.chrom, settings.exclude_sites_fn, ref);
  }

  // if(false){
  //   std::cerr << "\t-> TEMPORARILY:  REMOVING LOW COVERAGE CPGs from REF" << '\n';
  //   std::string fn ("temp/meth_low_data_to_be_removed.txt");
  //   filter_ref_db(settings.chrom, fn, ref);
  // } else {
  //   std::cerr << "\t-> TEMPORARILY: NOT REMOVING LOW COVERAGE CPGs from REF" << '\n';
  // }

  std::vector<tally_mat_freq> tmf;
  tmf.resize(rgs.size());
  std::vector<tally_mat> tm;
  tm.resize(rgs.size());
  std::cerr << "\t-> Init Count matrix" << '\n';
  for (size_t i=0; i<rgs.size(); i++){
    tm[i] = init_tallymat(settings.max_pos_to_end);
  }

  my_cov_rg cov_rg(rgs.size());

  const keeplist_map cpg_map = get_cpg_genomic_pos(ref, seq_len);
  const std::vector<int> cpg_bool = get_cpg_bool(ref, seq_len);
  std::cerr << "\t-> " << cpg_map.size() << " CpG's in chrom: " << settings.chrom << '\n';

  std::vector<per_site> data;
  data.resize(cpg_map.size());
  size_t nocpg_to_include = 1e6;
  per_site_nocpg nocpg_data(nocpg_to_include); 

  bam1_t *rd = bam_init1();
  int reads = 1e6;
  int counter = 0;
  alignment_data d;
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

    // this is 50% faster than analyzing every single read.
    startpos = rd->core.pos;
    endpos = bam_endpos(rd);

    // auto hit = std::find(cpg_bool.begin()+startpos, cpg_bool.begin()+endpos, 1);
    // if(hit!=cpg_bool.begin()+endpos){
    //   skipread=false;
    // }
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

    if (rd->core.l_qseq < settings.minreadlength || rd->core.qual < settings.minmapQ || (rd->core.flag & settings.flags_off) != 0 || rd->core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY)){
      trashed++;
      continue;
    }

   // when adding the rd to data, we should also add the RG as a vector. It could be integers referecing the all the RG found in the header.. 

   if(rg_split){
     rgptr = bam_aux_get(rd, "RG");
     if(rgptr==NULL){
       // in the case norg exists. dont know if i should trash or put in UNKWOWN
       rname = std::string(bam_get_qname(rd));
       std::cerr << "[No_found_RG] :: " << rname<< '\n';
       trashed++;
       continue;
     } else {
       rgname = std::string( (const char*)(rgptr+1)); // removes the Z by +1. all colons are gone already
       auto match = std::find(rgs.begin(), rgs.end(), rgname);
       if(match != rgs.end()){
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
   add_aligned_data(settings, d, cpg_map, ref, data, tm[rgname_idx], nocpg_data, cov_rg, rgname_idx, cycles[rgname_idx]);
  }
  time(&end_time_load_data);
  std::cerr << "\t-> Processed: " << counter  << ". Reads filtered: " << 
    trashed << ". Reads skipped (nocpg overlap): " << reads_skipped << 
    ". Loaded in " << difftime(end_time_load_data, start_time_load_data) << " seconds." << '\n';
  
  for (size_t i=0; i<rgs.size(); i++){
    std::cerr << "\t-> Total_Observations RG: "<< rgs[i] << " CpG: " << cov_rg.cpg[i] << " CpGCoverage: " << (double)cov_rg.cpg[i]/(double)cpg_map.size() <<  ". CnoCpG: " << cov_rg.nocpg[i] << '\n';
  }
  for (size_t i=0; i<rgs.size(); i++){
    std::cerr << "\t-> Dumping count file: " << settings.outbase << "." << rgs[i] <<  ".tallycounts" << '\n';
    dump_count_file(settings, tm[i], rgs[i]);
    tmf[i] = read_count_file(settings, rgs[i]);
  }

  std::vector<std::vector<double>> param_deam;
  param_deam.resize(rgs.size());

  for (size_t i=0; i<rgs.size(); i++){
    // this fucks up the order of the files. 
    // its better to just remove sites with missing data, when merging alle the RGS
    // filter_per_site_data(settings, tmf[i], data[i], rgs[i]);

    if(!settings.deamrates_filename.empty() && check_file_exists(settings.deamrates_filename)){
      std::cerr << "\t-> Loading deamination rates from " << settings.deamrates_filename << '\n';
      std::cerr << "\t-> Make sure that the file contains the same number of pos to include. Deammeth does not check that" << '\n';
      param_deam[i] = load_deamrates_f(settings);
    } else if (check_file_exists(settings.outbase+"."+rgs[i]+".deamrates")){
      std::cerr << "\t-> Loading deamination rates from " << settings.outbase+"."+rgs[i]+".deamrates" << '\n';
      param_deam[i] = load_deamrates(settings, rgs[i]);
    }else {
      std::cerr << "\t-> Starting Optim of deamination rates. RG: " << rgs[i] << '\n';
      run_deamrates_optim(settings, tmf[i], data, nocpg_data, rgs[i], i, cov_rg);
      std::cerr << "\t-> Dumping deamination parameters to " << settings.outbase+"."+rgs[i]+".deamrates" << '\n';
      param_deam[i] = load_deamrates(settings, rgs[i]);
    }
  }
  std::vector<pre_calc_per_site> mle_data;

  double deamin_methylated, deamin_unmethylated;
  double noM, M;
  std::cerr << "\t-> Precalculating many things for ML of F." << '\n';
  mle_data.reserve(data.size());
  for (const auto &site : data) {
    if(!site.depth){
      continue;
    }
    pre_calc_per_site d(site);
    for (size_t i = 0; i < site.depth; i++) {
      deamin_methylated = param_deam[site.rgs[i]][get_param_idx(settings.max_pos_to_end, METHSTATE, site.pos_to_end[i], site.prime[i])];
      deamin_unmethylated = param_deam[site.rgs[i]][get_param_idx(settings.max_pos_to_end, UNMETHSTATE, site.pos_to_end[i], site.prime[i])];
      noM = calc_prob_obs_base(site.seqerrors[i], site.base_compos[i], deamin_unmethylated);
      M = calc_prob_obs_base(site.seqerrors[i], site.base_compos[i], deamin_methylated);

      d.pre_noM.push_back(noM);
      d.pre_M.push_back(M);

      d.maperrors.push_back(site.maperrors[i]);

      d.prob_no_deamin_event_meth *= 1-deamin_methylated;
      d.prob_no_deamin_event_unmeth *= 1-deamin_unmethylated;
      d.summary[site.base_compos[i]]++;
    }
    mle_data.push_back(d);
  }
  data.clear();
  std::cerr << "\t-> Dumping MLE of F to " << settings.outbase << ".F" << '\n';
  run_mle(settings, mle_data);
  std::cerr << "\t-> Cleaning up." << '\n';
  mle_data.clear();
  free(ref);
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(rd);
  bam_hdr_destroy(header);
  sam_close(in);
  time(&end_time);
  double time_gone = difftime(end_time, start_time);
  size_t minutes = time_gone / 60.0;
  size_t seconds = (int)time_gone % 60;
  std::cerr << "\t-> Done in " << minutes << ":" << seconds << " M:S." << '\n';

  return 0;
}

int main(int argc, char *argv[]) {
  return parse_bam(argc, argv);

}
