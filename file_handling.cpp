#include <iostream> // stdout/stdin/stderr
#include <string>
#include <vector>
#include <sstream> //stringstream
#include <fstream>  // open file
#include <utility> // pair


#include "file_handling.hpp"

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

void filter_ref_sites(const std::string & selected_chrom, std::string & filename, char * ref){
   std::ifstream f (filename.c_str());
   checkfilehandle<std::ifstream>(f, filename);
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
  checkfilehandle<std::ifstream>(f, filename);
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




// funcs for args parsing
std::vector<std::string> parse_chrom_file(std::string & filename){
  std::vector<std::string> res;
  std::ifstream f (filename.c_str());
  checkfilehandle<std::ifstream>(f, filename);
  if(f.is_open()){
    std::string row, chrom;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> chrom;
      res.push_back(chrom);
      ss.clear();
    }
  }
  return res;
}


std::vector<std::pair<size_t, size_t>> parse_bed_file(std::string filename,  std::string &curr_chrom){
  std::vector<std::pair<size_t, size_t>> res;
  std::ifstream f (filename.c_str());
  checkfilehandle<std::ifstream>(f, filename);
  if(f.is_open()){
    std::string row;
    std::string chrom;
    size_t start,end ;
    while(getline(f, row)){
      if(row.empty()){
        std::cerr << "BED: "<< filename << " contains empty rows. EXITING" << '\n';
        exit(EXIT_FAILURE);
      }
      std::stringstream ss(row);

      ss >> chrom >> start >> end;
      if(start>=end){
        std::cerr << "BED: "<< filename << " is not normal: chr: " << chrom << " Start: " << start << " END: " << end << ". EXITING" << '\n';
        exit(EXIT_FAILURE);
      }
      if(chrom==curr_chrom){
        res.push_back(std::make_pair (start,end));
      }
    }
  }
  f.close();
  return(res);
}
