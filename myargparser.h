#include <string>
#include <limits>
#include <fstream>
struct general_settings{
  std::ofstream args_stream;
  std::string bam_fn, reference_fn, chrom;
  std::string exclude_sites_fn, all_options, outbase, deamrates_filename;
  std::string readgroups_f, bed_f, priors_str;
  size_t minmapQ, minbaseQ, max_pos_to_end;
  double M;
  int flags_off;
  size_t cycles, minreadlength;
  size_t windowsize, max_cpgs;
  general_settings(){
    minmapQ = 25;
    max_pos_to_end = 125;
    minbaseQ = 20;
    minreadlength=25;
    flags_off = 3844; // https://broadinstitute.github.io/picard/explain-flags.html
    M = 0.75; // ~80% of all CpG are methylated in humans
    cycles=std::numeric_limits<size_t>::max();
    windowsize=std::numeric_limits<size_t>::max();
    max_cpgs=std::numeric_limits<size_t>::max();
  }
};

void args_parser( int argc, char* argv[], general_settings & settings );
