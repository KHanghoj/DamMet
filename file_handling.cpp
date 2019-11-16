#include <iostream> // stdout/stdin/stderr
#include <string>
#include <vector>
#include <sstream> //stringstream
#include <fstream>  // open file
#include <string>

#include "file_handling.hpp"

bool check_file_exists(std::string filename){
  std::ifstream f(filename.c_str(), std::ios::in);
  return f.good();
}


template <typename T>
void checkfilehandle(T &fh, std::string filename){
  if (! fh.is_open()){
    std::cerr << "Couldnt open file: " << filename << " EXITING " << std::endl;
    exit(EXIT_FAILURE);
  }
}

void filter_ref_sites(const std::string & selected_chrom, std::string & filename, char * ref){
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
        // -1 as it has to be zero-based :)
        if(pos>0){
          ref[pos-1] = 'N';
        } else {
          std::cerr << "Something is not correct in sites file: " << filename
                    << " at region: " << chrom << " " << pos
                    << " EXITING!!!!! " << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      ss.clear();
    }
  }
}


void filter_ref_bed(const std::string & selected_chrom, std::string & filename, char * ref){
  std::ifstream f (filename.c_str());
  checkfilehandle(f, filename);
  if(f.is_open()){
    std::string row, chrom;
    size_t start, end;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> chrom >> start >> end;
      if (chrom == selected_chrom) {
        if (start >= end) {
          std::cerr << "Something is not correct in BED file: " << filename
                    << " at region: " << chrom << " " << start << " " << end
                    << " EXITING!!!!! " << std::endl;
          exit(EXIT_FAILURE);
        }
        for (size_t val = start; val < end; val++) {
          ref[val] = 'N';
        }
      }
      ss.clear();
    }
  }
  f.close();
}
