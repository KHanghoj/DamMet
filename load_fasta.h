#include "nuclconv.h"
#include <htslib/faidx.h>
#include <htslib/hts.h>

// clang++ load_fasta.cc  -DMYDEBUG ./htslib/libhts.a -I./htslib/ -lz -lm -lbz2 -llzma -o load_fasta && time ./load_fasta ~/devel/test_angsd/FASTA/humans/hs37d5.chrom20.fa

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
#ifdef MYDEBUG
  int j=100000, maxprint =3;
  std::cout << "Contig " << chrom << " is " << seq_len << " long. Printing nucleotides (0-index) " << j << " to " << j + maxprint<< '\n';
  std::cout << "Correspond to samtools faidx ./fasta.fa " <<  " "  << chrom << ":" << j+1 << "-" << j + maxprint+1<< '\n';
  for (int i=0; i <= maxprint; i++){
    std::cout << ref[i+j] << " " << refToInt[(int)(ref[i+j])] << '\n';
  }
#endif
  return ref;
}
