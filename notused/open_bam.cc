#include "load_fasta.h"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <cmath>

// clang++ open_bam.cc -std=c++11  -DMYDEBUG ./htslib/libhts.a -I./htslib/ -lz -lm -lbz2 -llzma -o open_bam && time ./open_bam ~/devel/test_angsd/FASTA/humans/hs37d5.chrom20.fa ~/devel/test_angsd/BAM/humans/MethSimWithOutDamageWithOutSeqErrors.chrom20.X40.bam 20

double phred_to_double(const int & phred){
  return std::pow(10, -(phred/10.0));
}

int parse_bam(int argc, char * argv[]) {
  std::string ref_filename (argv[1]);
  std::string bam_filename (argv[2]);
  std::string chrom (argv[3]);  // 20
  // std::string chrom ( "20" );
  if(!bam_filename.empty()) {
    //open BAM for reading
    samFile *in = sam_open(bam_filename.c_str(), "r");
    if(in == NULL) {
      std::cerr << "Unable to open BAM/SAM file: " << bam << '\n';
      exit(EXIT_FAILURE);
    }

    //Get the header
    bam_hdr_t *header = sam_hdr_read(in);
    if(header == NULL) {
      sam_close(in);
      std::cerr << "Unable to open BAM header: " << bam << '\n';
    }

    //Load the index
    hts_idx_t *idx = sam_index_load(in, bam_filename.c_str());
    if(idx == NULL) {
      std::cerr << "Unable to open BAM/SAM index. Make sure alignments are indexed." << '\n';
    }

    // go to the correct chromosome
    // https://github.com/gatoravi/bam-parser-tutorial/blob/master/parse_bam.cc
    hts_itr_t * iter  = sam_itr_querys(idx, header, chrom.c_str());

    // load ref of that chromosome
    faidx_t *fai = ref_init(ref_filename);
    char * ref = fetch_chrom(fai, chrom);

    //Initiate the alignment record
    bam1_t *aln = bam_init1();
    while(sam_itr_next(in, iter, aln) >= 0) {

      std::cout << "Read Chr: " << header->target_name[aln->core.tid];
      std::cout << "\tPos: " << aln->core.pos;
      std::vector<int> seqv, qualv;
      std::string seq, qual;
      std::vector<int> ref_read;
      uint8_t *quali = bam_get_qual(aln);
      uint8_t *seqi = bam_get_seq(aln);
      for (int i = 0; i < aln->core.l_qseq; i++) {
        ref_read.push_back(refToInt[(int)ref[i+aln->core.pos]]);
        seqv.push_back(seq_nt16_int[bam_seqi(seqi, i)]);
        seq += seq_nt16_str[bam_seqi(seqi, i)];
        qualv.push_back(quali[i]);
        qual += 33 + quali[i];
      }
      std::cout << "\tSeq: " << seq << "\tQual: " << qual << '\n';
      std::cout << "Seq:";
      for (const auto & val: seqv){
        std::cout << " " << val;
      }
      std::cout << "\nRef:";
      for (const auto & val: ref_read){
        std::cout << " " << val;
      }
      std::cout << "\nQual:";
      for (const auto & val: qualv){
        std::cout << " " << val;
      }
      std::cout << "\nQual:";
      for (const auto & val: qualv){
        std::cout << " " << phred_to_double(val);
      }
      std::cout << std::endl;
      break;
    }

    hts_itr_destroy(iter);
    hts_idx_destroy(idx);
    bam_destroy1(aln);
    bam_hdr_destroy(header);
    sam_close(in);
  }
  return 0;
}


int main(int argc, char* argv[]) {
  return parse_bam(argc, argv);
}
