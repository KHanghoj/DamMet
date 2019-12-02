#pragma once
#include <string.h>
#include <string>
#include <vector>
#include <limits>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h> // mkdir
#include <sstream> //stringstream

// #include "file_handling.hpp"

bool check_file_exists(std::string filename);
std::vector<std::string> parse_chrom_file(std::string & filename);

struct general_settings{
  std::ofstream args_stream;
  std::string buffer;
  std::string analysis;
  std::string bam_fn, reference_fn, chrom_temp;
  std::string exclude_sites_fn, exclude_bed_fn, outbase, deamrates_filename;
  std::string readgroups_f, bed_f, priors_str;
  std::vector<std::string> chrom;
  bool skip_empty_cpg=false;
  bool verbose=false;
  size_t minmapQ, minbaseQ, max_pos_to_end;
  double M;
  int flags_off;
  size_t cycles, minreadlength, minreadlength_deam;
  size_t windowsize, max_cpgs;
  general_settings(){
    minmapQ = 25;
    max_pos_to_end = 30;
    minbaseQ = 20;
    minreadlength=25;
    minreadlength_deam=std::numeric_limits<size_t>::max();
    flags_off = 3844; // https://broadinstitute.github.io/picard/explain-flags.html
    M = 0.75; // ~70-80% of all CpG are methylated in mammals
    cycles=std::numeric_limits<size_t>::max();
    windowsize=std::numeric_limits<size_t>::max();
    max_cpgs=std::numeric_limits<size_t>::max();
    buffer = ("");
  }
};

void args_parser( int argc, char* argv[], general_settings & settings );
void print_help();
bool check_estF_args(general_settings & settings);
