#include "myargparser.hpp"
#include "version.hpp"


void print_help(){
  std::cerr << "\nDamMet (" << DAMMET_VERSION << ") is a software aimed to estimate methylation maps using HTS sequencing "
    "data underlying ancient samples. The implemented model follows a two-steps procedure. "
    "The first step obtains a Maximum Likelihood Estimate (MLE) of position-specific deamination "
    "rates at both methylated and unmethylated cytosine residues. The second step makes use "
    "of these estimates to recover a MLE of local methylation levels in a user-defined window size." << std::endl;
  std::cerr << "Three args are required:\n\t->-b (bam)\n\t->-r (reference "
    "fasta)\n\t->-c (chromosome of interest)" << '\n';
  std::cerr << "OPTIONS:" << std::endl;
  std::cerr << "\t-> BED file (-B): " << std::endl;
  std::cerr << "\t-> minmapQ (-q): " << std::endl;
  std::cerr << "\t-> minbaseQ (-Q): " << std::endl;
  std::cerr << "\t-> MinReadLength (-L): " << std::endl;
  std::cerr << "\t-> MinReadLength_Deamrates (-l): " << std::endl;
  std::cerr << "\t-> Max_Pos_From_End (-P): " << std::endl;
  std::cerr << "\t-> Expected fraction of methylated CpGs (-M): " << std::endl;
  std::cerr << "\t-> Outbase (-O): " << std::endl;
  std::cerr << "\t-> readFlags (-F): " << std::endl;
  std::cerr << "\t-> Number of cycles (-C) (Only used if no RG file is NOT provided): " << std::endl;
  std::cerr << "\t-> Using Precalc deamination rates from (-D): " << std::endl;
  std::cerr << "\t-> Using readgroups from file (-R): " << std::endl;
  std::cerr << "\t-> Exclude sites (1-based) (-E): " << std::endl;
  std::cerr << "\t-> Exclude BED (-e): " << std::endl;
  std::cerr << "\t-> WindowSize (-W): " << std::endl;
  std::cerr << "\t-> Max CpGs per Window (-N): " << std::endl;
}

std::vector<std::string> split_comma_chrom(std::string & chrom){
  std::vector<std::string> res;
  std::stringstream ss( chrom );
  while( ss.good() ){
    std::string substr;
    getline( ss, substr, ',' );
    res.push_back( substr );
  }
  return res;
}

std::vector<std::string> parse_chrom_str(std::string & d){
  if (d.find(",") != std::string::npos){
    return split_comma_chrom(d);
  } else {
    std::vector<std::string> res;
    res.push_back(d);
    return(res);
  }
}
  
void args_parser(int argc, char *argv[], general_settings & settings) {
  if( (argc== 1) ||
      (argc== 2 && std::string(argv[1]) == "-h") ||
      (argc== 2 && std::string(argv[1]) == "-help") ||
      (argc== 2 && std::string(argv[1]) == "--help") ){
    print_help();
    exit(EXIT_SUCCESS);
  }

  std::vector<std::string> args;
  for(int i=1;i<argc;i++){
    args.push_back(std::string(argv[i]));
  }

  std::string ss("");
  for(auto i=args.begin();i!=args.end();i++){
    if(*i == "-bam" || *i == "-b" ){
      i++;
      if(check_file_exists(*i)){
        settings.bam_fn = *i;
      }else{
        std::cerr << "BAM FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
      ss += "\t-> BAM (-b): "+settings.bam_fn+'\n';
        
    }

    if((*i) == "-ref" || *i == "-r"){
      i++;
      if(check_file_exists(*i)){
        settings.reference_fn = *i;
      }else{
        std::cerr << "REFERENCE FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
      ss += "\t-> REF (-r): " + settings.reference_fn + '\n';
    }

    if((*i) == "-c" || *i == "-chrom"){
      i++;
      settings.chrom = parse_chrom_str(*i);
      ss += "\t-> Chromosome (-c): " + (*i) + '\n';
    }

    if((*i) == "-cf" || *i == "--chromfile"){
      i++;
      if(check_file_exists(*i)){
         settings.chrom = parse_chrom_file(*i);
      }else{
        std::cerr << "CHROM FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
      ss += "\t-> Chromosome (-cf): " + (*i) + '\n';
    }
 

    if((*i) == "-O"){
      i++;
      settings.outbase = *i;
    }


    if((*i) == "-M"){
      i++;
      settings.M = std::stof(*i);
      if(settings.M < 0 || settings.M > 1){
        std::cerr << "Fraction of methylated CpGs on the chromosome must be a float between 0 and 1" << '\n';
        exit(EXIT_FAILURE);
      }
    }

    if((*i) == "-W"){
      i++;
      settings.windowsize = std::stol(*i);
      ss += "\t-> WindowSize (-W): " + std::to_string(settings.windowsize) + '\n';
    }

    if((*i) == "-N"){
      i++;
      settings.max_cpgs = std::stoi(*i);
      ss += "\t-> Max CpGs per Window (-N): " + std::to_string(settings.max_cpgs) + '\n';
    }

    if((*i) == "-B"){
      i++;
      if(check_file_exists(*i)){
        settings.bed_f = *i;
      }else{
        std::cerr << "BED FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
      
      ss += "\t-> BED file (-B): " + settings.bed_f + '\n';      
    }
   
    if((*i) == "-q"){
      i++;
      settings.minmapQ = std::stoi(*i);
      if(settings.minmapQ < 1){
        settings.minmapQ=1;
      }

    }

    if((*i) == "-Q"){
      i++;
      settings.minbaseQ = std::stoi(*i);
      if(settings.minbaseQ < 1){
        settings.minbaseQ=1;
      }

    }

    if((*i) == "-F"){
      i++;
      settings.flags_off = std::stoi(*i);
    }

    if((*i) == "-C"){
      i++;
      settings.cycles = std::stoi(*i);
    }

    if((*i) == "-P"){
      i++;
      // 5 -> 0,1,2,3,4. so 5-1
      settings.max_pos_to_end = std::stoi(*i)-1;
    }
    
    if((*i) == "-E" ){
      i++;
      if(check_file_exists(*i)){
        settings.exclude_sites_fn = *i;
      }else{
        std::cerr << "SITES FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
    }

    if((*i) == "-e" ){
      i++;
      if(check_file_exists(*i)){
        settings.exclude_bed_fn = *i;
      }else{
        std::cerr << "BED FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
    }
    
    if((*i) == "-D" ){
      i++;
      if(check_file_exists(*i)){
        settings.deamrates_filename = *i;
      }else{
        std::cerr << "DEAMRATES FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
      
    }

    if((*i) == "-h" ){
      i++;
      settings.priors_str = *i;
    }

    if((*i) == "-L" ){
      i++;
      settings.minreadlength = std::stoi(*i);
    }

    if((*i) == "-l" ){
      i++;
      settings.minreadlength_deam = std::stoi(*i);
    }

    if((*i) == "-R" ){
      i++;
      if(check_file_exists(*i)){
        settings.readgroups_f = *i;
      }else{
        std::cerr << "readgroups_f FILE DOES NOT EXISTS" << std::endl;
        exit(EXIT_FAILURE);        
      }
    }
    
  }


  if (settings.bam_fn.empty() || settings.reference_fn.empty() ||
      settings.chrom.empty()) {
    print_help();
    exit(EXIT_FAILURE);
  }

  if(settings.outbase.empty()){
    settings.outbase="dammet_res";
    mkdir(settings.outbase.c_str(),0777);
    settings.outbase += "/dammet";
    std::cerr << "\t-> output prefix (-O) is set to: " << settings.outbase << '\n';
  }

  
  if(settings.max_cpgs==std::numeric_limits<size_t>::max() && settings.windowsize==std::numeric_limits<size_t>::max() && settings.bed_f.empty()){
    print_help();
    std::cerr << "Must specify either -N (max CpGs per window) AND/OR -W (max windowsize) OR -B bedfile" << '\n';
    std::cerr << "EXITING...." << '\n';
    exit(EXIT_FAILURE);
  }

  if(settings.minreadlength_deam==std::numeric_limits<size_t>::max()){
    settings.minreadlength_deam=settings.minreadlength;
  }


  ss += "\t-> minmapQ (-q): " + std::to_string(settings.minmapQ) + '\n';
  ss += "\t-> minbaseQ (-Q): " + std::to_string(settings.minbaseQ) + '\n';
  ss += "\t-> MinReadLength (-L): " + std::to_string(settings.minreadlength) + '\n';
  ss += "\t-> MinReadLength_deamrates (-l): " + std::to_string(settings.minreadlength_deam) + '\n';
  ss += "\t-> Max_Pos_From_End (-P): " + std::to_string(settings.max_pos_to_end) + '\n';
  ss += "\t-> Expected fraction of methylated CpGs (-M): " + std::to_string(settings.M) + '\n';
  ss += "\t-> Outbase (-O): " + settings.outbase + '\n';
  ss += "\t-> readFlags (-F): " + std::to_string(settings.flags_off) + '\n';
  ss += "\t-> Number of cycles (-C) (Only used if no RG file is NOT provided): " + std::to_string(settings.cycles) + '\n';
  ss += "\t-> Using Precalc deamination rates from (-D): " + settings.deamrates_filename + '\n';
  ss += "\t-> Using readgroups from file (-R): " + settings.readgroups_f + '\n';
  ss += "\t-> Exclude sites (1-based) (-E): " + settings.exclude_sites_fn + '\n';
  ss += "\t-> Exclude BED (-e): " + settings.exclude_bed_fn + '\n';  
  
  settings.all_options = ss;
  // std::cerr << '\n' << ss;
}
