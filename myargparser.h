#include <string>
#include <limits>
struct general_settings{
  std::string bam_fn, reference_fn, chrom;
  std::string exclude_sites_fn, outbase, all_options, deamrates_filename;
  std::string deamrates_f_rg;
  size_t minmapQ, minbaseQ, max_pos_to_end;
  double het_rate;
  double M;
  int flags_off;
  size_t cycles, minreadlength;
  size_t windowsize, max_cpgs;
  general_settings(){
    het_rate = 0.001;
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

general_settings args_parser( int argc, char* argv[] );
