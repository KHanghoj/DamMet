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
