#include "myargparser.hpp"
#include <cstdlib>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h> // mkdir



void args_parser(int argc, char *argv[], general_settings & settings) {
  int c;
  // : means that an options is expected

  while ((c = getopt(argc, argv, "b:r:c:q:Q:F:E:e:C:P:D:h:O:M:W:N:L:R:B:l:")) != -1) {
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
    case 'e':
      if (optarg)
        settings.exclude_bed_fn = std::string(optarg);
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
        settings.priors_str = std::string(optarg);
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
    case 'l':
      if (optarg){
        // windowsize
        settings.minreadlength_deam = atoi(optarg);
      }
      break;
    case 'R':
      if (optarg){
        // one RG per line file
        settings.readgroups_f = std::string(optarg);
      }
      break;
    case 'B':
      if (optarg){
        // one RG per line file
        settings.bed_f = std::string(optarg);
      }
      break;
    case '?':
      std::cerr << "Unknown argument" << '\n';
      exit(EXIT_FAILURE);
      break;
    }
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

  if(settings.minmapQ < 1){
    settings.minmapQ=1;
  }
  if(settings.minbaseQ < 1){
    settings.minbaseQ=1;
  }

  if(settings.minreadlength_deam==std::numeric_limits<size_t>::max()){
    settings.minreadlength_deam=settings.minreadlength;
  }

  std::string ss("");
  ss += "\t-> BAM (-b): "+settings.bam_fn+'\n';
  ss += "\t-> REF (-r): " + settings.reference_fn + '\n';
  ss += "\t-> Chromosome (-c): " + settings.chrom + '\n';

  ss += "\t-> BED file (-B): " + settings.bed_f + '\n';
  ss += "\t-> minmapQ (-q): " + std::to_string(settings.minmapQ) + '\n';
  ss += "\t-> minbaseQ (-Q): " + std::to_string(settings.minbaseQ) + '\n';
  ss += "\t-> MinReadLength (-L): " + std::to_string(settings.minreadlength) + '\n';
  ss += "\t-> MinReadLength_deamrates (-l): " + std::to_string(settings.minreadlength_deam) + '\n';
  ss += "\t-> Max_Pos_From_End (-P): " + std::to_string(settings.max_pos_to_end+1) + '\n';
  ss += "\t-> Expected fraction of methylated CpGs (-M): " + std::to_string(settings.M) + '\n';
  ss += "\t-> Outbase (-O): " + settings.outbase + '\n';
  ss += "\t-> readFlags (-F): " + std::to_string(settings.flags_off) + '\n';
  ss += "\t-> Number of cycles (-C) (Only used if no RG file is NOT provided): " + std::to_string(settings.cycles) + '\n';
  ss += "\t-> Using Precalc deamination rates from (-D): " + settings.deamrates_filename + '\n';
  ss += "\t-> Using readgroups from file (-R): " + settings.readgroups_f + '\n';
  ss += "\t-> Exclude sites (1-based) (-E): " + settings.exclude_sites_fn + '\n';
  ss += "\t-> Exclude BED (-e): " + settings.exclude_bed_fn + '\n';  
  ss += "\t-> WindowSize (-W): " + std::to_string(settings.windowsize) + '\n';
  ss += "\t-> Max CpGs per Window (-N): " + std::to_string(settings.max_cpgs) + '\n';
  
  settings.all_options = ss;
  // std::cerr << '\n' << ss;
}
