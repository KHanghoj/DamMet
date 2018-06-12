#include "myargparser.h"
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h> // mkdir

general_settings args_parser(int argc, char *argv[]) {

  general_settings settings;
  int c;
  // : means that an options is expected

  while ((c = getopt(argc, argv, "b:r:c:q:Q:F:E:C:P:D:h:O:M:W:N:L:R:")) != -1) {
    switch (c) {
    case 'b':
      if (optarg)
        settings.bam_fn = std::string(optarg);
      break;
    case 'r':
      if (optarg)
        settings.reference_fn = std::string(optarg);
      break;
    case 'c':
      if (optarg)
        settings.chrom = std::string(optarg);
      break;
    // case 'T':
    //   if (optarg){
    //     settings.tally_on = atoi(optarg);
    //   }
    //   break;
    case 'q':
      if (optarg)
        settings.minmapQ = atoi(optarg);
      break;
    case 'Q':
      if (optarg)
        settings.minbaseQ = atoi(optarg);
      break;
    case 'F':
      if (optarg)
	// settings.flags_off = 0;
        // settings.flags_off = strtol(optarg, 0, 0); // = atoi(optarg);
	settings.flags_off = atoi(optarg);
      break;
    case 'E':
      if (optarg)
        settings.exclude_sites_fn = std::string(optarg);
      break;
    case 'C':
      if (optarg)
        settings.cycles = atoi(optarg);
      break;
    case 'P':
      // 5 -> 0,1,2,3,4. so 5-1
      if (optarg)
        settings.max_pos_to_end = atoi(optarg)-1;
      break;
    case 'D':
      if (optarg){
        settings.deamrates_filename = std::string(optarg);
      }
      break;
    case 'h':
      if (optarg){
        settings.het_rate = atof(optarg);
      }
      break;
    case 'O':
      if (optarg){
        settings.outbase = std::string(optarg);
      }
      break;
    case 'M':
      if (optarg){
        // expected fraction of methylated CpG's
        settings.M = atof(optarg);
      }
      break;
    case 'W':
      if (optarg){
        // windowsize
        settings.windowsize = atol(optarg);
      }
      break;
    case 'N':
      if (optarg){
        // windowsize
        settings.max_cpgs = atoi(optarg);
      }
      break;
    case 'L':
      if (optarg){
        // windowsize
        settings.minreadlength = atoi(optarg);
      }
      break;
    case 'R':
      if (optarg){
        // one RG per line file
        settings.deamrates_f_rg = std::string(optarg);
      }
      break;
    case '?':
      std::cerr << "FUCKED BIGTIME" << '\n';
      exit(EXIT_FAILURE);
      break;
    }
  }

  if (settings.bam_fn.empty() || settings.reference_fn.empty() ||
      settings.chrom.empty()) {
    std::cerr << "Three args are required:\n\t->-b (bam)\n\t->-r (reference "
                 "fasta)\n\t->-c (chromosome of interest)\n\t->-C (read/dump "
                 "count table)"
              << '\n';
    exit(EXIT_FAILURE);
  }

  if(settings.M < 0 || settings.M > 1){
    std::cerr << "Fraction of methylated CpGs on the chromosome must be between 0 and 1" << '\n';
    exit(EXIT_FAILURE);
  }

  if(settings.outbase.empty()){
    settings.outbase="dammet_res";
    mkdir(settings.outbase.c_str(),0777);
    settings.outbase += "/"+settings.chrom;
  }

  if(settings.max_cpgs==std::numeric_limits<size_t>::max() && settings.windowsize==std::numeric_limits<size_t>::max() ){
    std::cerr << "Must specify either -N (max CpGs per window) AND/OR -W (max windowsize)" << '\n';
    std::cerr << "EXITING...." << '\n';
    exit(EXIT_FAILURE);
  }

  if(settings.het_rate > 1){
    std::cerr << "Heterozygosity prior must be between 0 and 1. Can be set to -1 for uniform prior" << '\n';
    exit(EXIT_FAILURE);
  }

  if(settings.het_rate < 0){
    std::cerr << "using a uniform prior for heterozygosity estimates" << '\n';
    settings.het_rate=-1;
  }

  if(settings.minmapQ < 1){
    settings.minmapQ=1;
  }
  if(settings.minbaseQ < 1){
    settings.minbaseQ=1;
  }

  std::string ss("");
  ss += "\t-> BAM (-b): "+settings.bam_fn+'\n';
  ss += "\t-> REF (-r): " + settings.reference_fn + '\n';
  ss += "\t-> Chromosome (-c): " + settings.chrom + '\n';
  ss += "\t-> minmapQ (-q): " + std::to_string(settings.minmapQ) + '\n';
  ss += "\t-> minbaseQ (-Q): " + std::to_string(settings.minbaseQ) + '\n';
  ss += "\t-> MinReadLength (-L): " + std::to_string(settings.minreadlength) + '\n';
  ss += "\t-> Max_Pos_From_End (-P): " + std::to_string(settings.max_pos_to_end+1) + '\n';
  ss += "\t-> Expected fraction of methylated CpGs (-M): " + std::to_string(settings.M) + '\n';
  ss += "\t-> Outbase (-O): " + settings.outbase + '\n';
  ss += "\t-> readFlags (-F): " + std::to_string(settings.flags_off) + '\n';
  ss += "\t-> Number of cycles (-C) (Only used if no RG file is NOT provided): " + std::to_string(settings.cycles) + '\n';
  if(!settings.deamrates_filename.empty()){
    ss += "\t-> Using Precalc deamination rates from (-D): " + settings.deamrates_filename + '\n';
  }
  ss += "\t-> Using readgroups from file (-R): " + settings.deamrates_f_rg + '\n';
  ss += "\t-> Exclude sites (-E): " + settings.exclude_sites_fn + '\n';
  
  ss += "\t-> WindowSize (-W): " + std::to_string(settings.windowsize) + '\n';
  ss += "\t-> Max CpGs per Window (-N): " + std::to_string(settings.max_cpgs) + '\n';
  
  settings.all_options = ss;
  std::cerr << '\n' << ss;
  return settings;
}
