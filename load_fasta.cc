#include "load_fasta.h"

// extern const int seq_nt16_int[] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };
// clang++ load_fasta.cc  -DMYDEBUG ./htslib/libhts.a -I./htslib/ -lz -lm -lbz2 -llzma -o load_fasta && time ./load_fasta test_angsd/FASTA/humans/hs37d5.chrom20.fa

int main(int argc, char *argv[])
{
  int optind = 1;
  std::string filename (argv[optind]);
  // char * chrom = argv[2];
  std::string chrom ("20");
  faidx_t *fai = ref_init(filename);
  char * ref = fetch_chrom(fai, chrom);
  // ref = fai_fetch(fai, h->target_name[tid], &seq_len);
  free(ref);
  fai_destroy(fai);
  return 0;
}
