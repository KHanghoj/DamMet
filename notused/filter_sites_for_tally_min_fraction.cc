#include <fstream>
#include <sstream> //stringstream
#include "load_fasta.h"
#include <htslib/hts.h>
#include <htslib/sam.h>

size_t NUCLEOTIDES =4 ;

const std::string DELIMITERS = " \t" ;

inline void getsiteinfo(const std::string &row, size_t &lastOffset, std::string &result){
  int offset = row.find_first_of(DELIMITERS, lastOffset);
  result = row.substr(lastOffset, offset - lastOffset);
  lastOffset = offset+1;
}

inline void get_next(const std::string &row, size_t &lastOffset, size_t &result){
  int offset = row.find_first_of(DELIMITERS, lastOffset);
  result = atoi(row.substr(lastOffset, offset - lastOffset).c_str());
  lastOffset = offset+1;
}

inline void get_next(const std::string &row, size_t &lastOffset, double &result){
  int offset = row.find_first_of(DELIMITERS, lastOffset);
  result = atof(row.substr(lastOffset, offset - lastOffset).c_str());
  lastOffset = offset+1;
}


int main(int argc, char *argv[])
{

  if (argc!=5 && argc != 6){
    std::cerr << "[REF] [CHROM] [MAX_ALLELE_COUNT] [OUTFILE]" << '\n';
    std::cerr << "[REF] [CHROM] [MAX_ALLELE_COUNT] [OUTFILE] [MAX_ALLOWED_FRACTION_ALT_ALLELE]" << '\n';
    exit(EXIT_FAILURE);
  }
  std::string ref_filename (argv[1]);
  std::string chrom (argv[2]);  // 20
  size_t max_alt_allele_count = atoi(argv[3]);
  std::string outfile (argv[4]);
  std::string outerrormat (outfile+".errormat");
  double max_allow_fraction_alt;
  if (argc==6){
    max_allow_fraction_alt = atof(argv[5]);
  }  else {
    max_allow_fraction_alt = 0.0;
  }
  std::cerr << max_allow_fraction_alt << '\n';

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
  size_t pos, A,C,G,T;
  double tot;
  std::stringstream ss;
  getline(std::cin, row); // header
  std::vector< std::vector< size_t > > ref_to_error (NUCLEOTIDES);
  for (auto & a: ref_to_error){
    a.resize(NUCLEOTIDES);
  }
  bool hetero_site = true;

  size_t offset;
  while(getline(std::cin, row)){
    offset=0;
    getsiteinfo(row, offset, curr_chrom);
    get_next(row, offset, pos);
    get_next(row, offset, tot);
    get_next(row, offset, A);
    get_next(row, offset, C);
    get_next(row, offset, G);
    get_next(row, offset, T);

    // ss.str(row);
    // ss >> curr_chrom >> pos >> tot >> A >> C >> G >> T;
    // ss.clear();
    size_t ref_base = refToInt[(int)ref[pos-1]];
    hetero_site = true;
    if (ref_base == 4) {
      continue;
    }

    if((ref_base == 0 && C+G+T<=max_alt_allele_count) || (ref_base == 0 && (double)(C+G+T) / tot <= max_allow_fraction_alt)){
      hetero_site = false;
    } else if((ref_base == 1 && A+G+T<=max_alt_allele_count) || (ref_base == 1 && (double)(A+G+T) / tot <= max_allow_fraction_alt)){
      hetero_site = false;
    } else if((ref_base == 2 && A+C+T<=max_alt_allele_count) || (ref_base == 2 && (double)(A+C+T) / tot <= max_allow_fraction_alt)){
      hetero_site = false;
    } else if((ref_base == 3 && A+C+G<=max_alt_allele_count) || (ref_base == 3 && (double)(A+C+G) / tot <= max_allow_fraction_alt)){
      hetero_site = false;
    }

    // if(ref_base == 0 && C+G+T<=max_alt_allele_count){
    //   hetero_site = false;
    // } else if(ref_base == 1 && A+G+T<=max_alt_allele_count){
    //   hetero_site = false;
    // } else if(ref_base == 2 && A+C+T<=max_alt_allele_count){
    //   hetero_site = false;
    // } else if(ref_base == 3 && A+C+G<=max_alt_allele_count){
    //   hetero_site = false;
    // }

    if (! hetero_site){
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
