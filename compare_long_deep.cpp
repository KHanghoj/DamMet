#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib> // std::div
#include <memory>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream> //stringstream



const size_t PRIMES=2;
const size_t STRANDS=2;
const size_t NUCLEOTIDES=4;

const size_t READPOS_MULT = PRIMES*STRANDS*std::pow(NUCLEOTIDES,4);
const size_t PRIME_MULT = STRANDS*std::pow(NUCLEOTIDES,4);
const size_t STRAND_MULT = std::pow(NUCLEOTIDES,4);
const size_t B1_MULT = std::pow(NUCLEOTIDES,3);
const size_t B2_MULT = std::pow(NUCLEOTIDES,2);
const size_t B3_MULT = NUCLEOTIDES;

using tally_mat = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<size_t>>>>>>>;

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
  //thematrix[1][0][0][1][2][3][2] = 100;
  return thematrix;
}

std::vector<size_t> long_vector(const int & readpos){
  int n = 1;
  // +1 to include the the array for alle obs within the reads
  n *= readpos+1;
  n *= PRIMES * STRANDS;
  n *= std::pow(NUCLEOTIDES,4);
  return std::vector<size_t>(n, 0);
}

size_t get_param_idx(const size_t & max_pos_to_end, const size_t & meth, const size_t & pos_to_end, const size_t & prime){
  return (meth * PRIMES * ( max_pos_to_end+2 ) + pos_to_end * PRIMES  + prime);
}

size_t get_idx(const size_t & readpos,
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

template <typename T>
void checkfilehandle(T &fh, std::string filename){
  if (! fh.is_open()){
    std::cerr << "Couldnt open file: " << filename << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
}


tally_mat read_file_deep(const std::string & fname, const int & n_readpos){
  std::ifstream f(fname.c_str());
  tally_mat temp = init_tallymat(n_readpos);
  if(f.is_open()){
    std::string row;
    size_t p, pr, strand, r1, r2, s1, s2, val;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> p >> pr >> strand >> r1 >> r2 >> s1 >> s2 >> val;
      temp[p][pr][strand][r1][r2][s1][s2] = val;
      ss.clear();
      }
    }
  f.close();
  return temp;
}

std::vector<size_t> read_file_long(const std::string & fname, const int & n_readpos){
  std::ifstream f(fname.c_str());
  std::vector<size_t> temp =  long_vector(n_readpos);;
  if(f.is_open()){
    std::string row;
    size_t p, pr, strand, r1, r2, s1, s2, val;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> p >> pr >> strand >> r1 >> r2 >> s1 >> s2 >> val;
      size_t idx = get_idx(p, pr, strand, r1, r2, s1, s2);
      temp[idx] = val;
      ss.clear();
      }
    }
  f.close();
  return temp;
}


void dump_count_file_deep(tally_mat & tm, const size_t & readpos){
  std::string filename = "out.deep.txt";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  if (f.is_open()) {
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

std::vector<int> idx_to_params(const int & idx){
  std::div_t p = std::div(idx, READPOS_MULT);
  std::div_t pr = std::div(p.rem, PRIME_MULT);
  std::div_t strand = std::div(pr.rem, STRAND_MULT);
  std::div_t B1 = std::div(strand.rem, B1_MULT);
  std::div_t B2 = std::div(B1.rem, B2_MULT);
  std::div_t B3 = std::div(B2.rem, B3_MULT);
  return std::vector<int> {p.quot, pr.quot, strand.quot, B1.quot, B2.quot, B3.quot, B3.rem};
}

void dump_count_file_long2(std::vector<size_t> & tm){
  std::string filename = "out.long2.txt";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  if (f.is_open()) {
    for(size_t n=0; n<tm.size(); n++){
      if(tm[n]){
        std::vector<int> info = idx_to_params(n);
         f << info[0] << " " << info[1] << " " << info[2] << " " << info[3] << " "
                      << info[4] << " " << info[5] << " " << info[6] << " "
                      << tm[n] << '\n';
      }
    }
    f.flush();
  }
  f.close();
}
void dump_count_file_long(std::vector<size_t> & tm, const size_t & readpos){
  std::string filename = "out.long.txt";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);
  if (f.is_open()) {
    for (size_t p = 0; p <= readpos; p++) {
      for (size_t pr = 0; pr < PRIMES; pr++) {
        for (size_t strand = 0; strand < STRANDS; strand++) {
          for (size_t r1 = 0; r1 < NUCLEOTIDES; r1++) {
            for (size_t r2 = 0; r2 < NUCLEOTIDES; r2++) {
              for (size_t s1 = 0; s1 < NUCLEOTIDES; s1++) {
                for (size_t s2 = 0; s2 < NUCLEOTIDES; s2++) {
                  size_t idx = get_idx(p, pr, strand, r1, r2, s1, s2);
                  if (tm[idx]) {
                    f << p << " " << pr << " " << strand << " " << r1 << " "
                      << r2 << " " << s1 << " " << s2 << " "
                      << tm[idx] << '\n';
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


int main(){
  int n_readpos = 2;
  std::vector<size_t> v = long_vector(n_readpos);
  std::cout << v.size() << '\n';


  size_t readpos = 1;
  size_t prime = 0;
  size_t strand = 0;
  size_t b1 = 1 ;
  size_t b2 = 2 ;
  size_t b3 = 3 ;
  size_t b4 = 2 ;

  // size_t idx = get_idx(readpos, prime, strand, b1, b2, b3, b4);
  // std::cerr << idx << '\n';
  // v[idx] = 100;

  tally_mat old_vec = read_file_deep("test.txt", 2);
  std::vector<size_t> new_vec = read_file_long("test.txt", 2);

  dump_count_file_deep(old_vec, 2);
  dump_count_file_long(new_vec, 2);
  dump_count_file_long2(new_vec);

  std::cerr << READPOS_MULT << '\n';
  std::cerr << PRIME_MULT << '\n';
  std::cerr << STRAND_MULT << '\n';
  std::cerr << B1_MULT << '\n';
  std::cerr << B2_MULT << '\n';
  std::cerr << B3_MULT << '\n';
  return 1;

}
