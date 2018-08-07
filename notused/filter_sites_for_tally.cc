#include <fstream>
#include <sstream> //stringstream
#include "load_fasta.h"
#include <htslib/hts.h>
#include <htslib/sam.h>

size_t NUCLEOTIDES =4 ;

int main(int argc, char *argv[])
{

  if (argc!=5){
    std::cerr << "[REF] [CHROM] [MAX_ALLELE_COUNT] [OUTFILE]" << '\n';
    exit(EXIT_FAILURE);
  }
  std::string ref_filename (argv[1]);
  std::string chrom (argv[2]);  // 20
  size_t max_alt_allele_count = atoi(argv[3]);
  std::string outfile (argv[4]);
  std::string outerrormat (outfile+".errormat");

  std::ofstream f(outfile.c_str());
  if(! f.is_open()){
    std::cerr << "Couldnt open file: " << outfile << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
  std::ofstream ferror(outerrormat.c_str());
  if(! ferror.is_open()){
    std::cerr << "Couldnt open file: " << outerrormat << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }


  faidx_t *fai = ref_init(ref_filename);
  char * ref = fetch_chrom(fai, chrom);
  // size_t seq_len = faidx_seq_len(fai, chrom.c_str());
  std::string row, curr_chrom;
  size_t pos, tot, A,C,G,T;
  std::stringstream ss;
  getline(std::cin, row); // header
  std::vector< std::vector< size_t > > ref_to_error (NUCLEOTIDES);
  for (auto & a: ref_to_error){
    a.resize(NUCLEOTIDES);
  }
  bool skip_site = false;
  while(getline(std::cin, row)){
    ss.str(row);
    ss >> curr_chrom >> pos >> tot >> A >> C >> G >> T;
    ss.clear();
    size_t ref_base = refToInt[(int)ref[pos-1]];
    skip_site = false;
    if (ref_base == 4) {
      continue;
    }

    if(ref_base == 0 && C+G+T<=max_alt_allele_count){
      skip_site = true;
    } else if(ref_base == 1 && A+G+T<=max_alt_allele_count){
      skip_site = true;
    } else if(ref_base == 2 && A+C+T<=max_alt_allele_count){
      skip_site = true;
    } else if(ref_base == 3 && A+C+G<=max_alt_allele_count){
      skip_site = true;
    }

    if (skip_site){
      ref_to_error[ref_base][0]+=A;
      ref_to_error[ref_base][1]+=C;
      ref_to_error[ref_base][2]+=G;
      ref_to_error[ref_base][3]+=T;
    } else {
      f << chrom << " " << pos << " " << intToRef[ref_base] << " " << tot  << " "  << A << " "  << C << " "  << G << " "  << T << '\n';
    }
  }

  for (size_t i=0; i<NUCLEOTIDES; i++){
    size_t total=0;
    for (size_t j=0; j<NUCLEOTIDES; j++){
      total+=ref_to_error[i][j];
    }
    for (size_t j=0; j<NUCLEOTIDES; j++){
      ferror << i << " " << j << " " << ref_to_error[i][j] << " " << ref_to_error[i][j]/(double)total <<  '\n';
    }
  }


  return 0;
}
