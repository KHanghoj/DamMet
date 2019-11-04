#include <iostream> // stdout/stdin/stderr
#include <string>
#include <vector>
#include "htslib/faidx.h"

faidx_t * ref_init(std::string & filename){
  faidx_t * fai = fai_load(filename.c_str());
  if ( !fai ) {
    std::cerr << "Could not load fai index of: " << filename << '\n';
    exit(EXIT_FAILURE);
  }
  return fai;
}

char * fetch_chrom(const faidx_t *fai, std::string & chrom)
{
  int seq_len = 0;
  char * ref = fai_fetch(fai, chrom.c_str(), &seq_len);
  if ( seq_len < 0 ) {
    std::cerr << "Failed to fetch sequence in: " <<  chrom << '\n';
    exit(EXIT_FAILURE);
  }
  return ref;
}
