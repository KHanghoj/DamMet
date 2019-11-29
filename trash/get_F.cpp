struct F_void {
  general_settings * settings;
  std::vector<pre_calc_per_site> * data;
  per_mle_run * mle_data;
  size_t iteration;
  F_void(general_settings * s, std::vector<pre_calc_per_site> * d, per_mle_run * m){
    settings = s;
    data = d;
    mle_data = m;
    iteration = 0;
  }
};



int est_dam_and_F(general_settings & settings) {

  // update priors if provided
  if(!settings.priors_str.empty()){
    update_dinucl_priors(settings);
  }

  std::string stream_filename (settings.outbase+".args");
  settings.args_stream.open(stream_filename.c_str());
  checkfilehandle(settings.args_stream, stream_filename);
  settings.args_stream << settings.all_options;


  // write priors to log
  std::cerr << "\t-> PRIORS (-h): " << std::exp(LOG_PRIORS[0]);
  settings.args_stream  << "\t-> PRIORS (-h): " << std::exp(LOG_PRIORS[0]);
  for (auto it=LOG_PRIORS.begin()+1; it!=LOG_PRIORS.end();it++){
    std::cerr  << "," << std::exp(*it);
    settings.args_stream  << "," << std::exp(*it);
  }
  std::cerr << '\n';
  settings.args_stream << '\n';


  std::vector<std::pair<size_t, size_t>> bed_coord;
  // parsing bedfile if provided
  if(!settings.bed_f.empty()){
    parse_bed_file(settings, bed_coord);
    std::cerr << "\t-> Analyzing: " << bed_coord.size() << " BED regions on chrom: " << settings.chrom << '\n';
    settings.args_stream << "\t-> Analyzing: " << bed_coord.size() << " BED regions on chrom: " << settings.chrom << '\n';

    if(bed_coord.size()==0){
      std::cerr << '\n' << "EXITING. NOT BED COORDS ON CHROM" << '\n';
      exit(EXIT_SUCCESS);
    }

  }

  // open BAM for reading
  samFile *in = sam_open(settings.bam_fn.c_str(), "r");
  if (in == NULL) {
    std::cerr << "Unable to open BAM/SAM file: " << bam << '\n';
    exit(EXIT_FAILURE);
  }

  // Get the header
  bam_hdr_t *header = sam_hdr_read(in);
  if (header == NULL) {
    sam_close(in);
    std::cerr << "Unable to open BAM header: " << bam << '\n';
    exit(EXIT_FAILURE);
  }

  // Load the index
  hts_idx_t *idx = sam_index_load(in, settings.bam_fn.c_str());

  if (idx == NULL) {
    std::cerr
        << "Unable to open BAM/SAM index. Make sure alignments are indexed."
        << '\n';
    exit(EXIT_FAILURE);
  }

  // load read groups
  bool rg_split;
  std::vector<std::string> rgs;
  std::vector<size_t> cycles;
  if((!settings.readgroups_f.empty()) && check_file_exists(settings.readgroups_f)){
    rg_split = true;
    std::ifstream f (settings.readgroups_f.c_str());
    checkfilehandle(f, settings.readgroups_f);
    std::string row, rg, cycle_string;
    size_t cycle;
    std::stringstream ss;
    while(getline(f, row)){
      ss.str(row);
      ss >> rg >> cycle_string;
      if(cycle_string.empty()){
        cycle=std::numeric_limits<size_t>::max();
        std::cerr << "\t-> RG: " << rg  << ". setting sequencing cycles: " << cycle << '\n';
        settings.args_stream << "\t-> RG: " << rg  << ". setting sequencing cycles: " << cycle << '\n';
      } else {
        cycle=std::stoi(cycle_string);
        std::cerr << "\t-> RG: " << rg  << ". sequencing cycles: " << cycle << '\n';
        settings.args_stream << "\t-> RG: " << rg  << ". sequencing cycles: " << cycle << '\n';
      }
      rgs.push_back(rg);
      cycles.push_back(cycle);
      ss.clear();
      rg.clear();
      cycle_string.clear();
    }
    f.close();
  } else {
    rg_split = false;
    rgs.push_back(ALL_RG);
    if(settings.cycles!=std::numeric_limits<size_t>::max()){
      cycles.push_back(settings.cycles);
    } else {
      cycles.push_back(std::numeric_limits<size_t>::max());
    }
    std::cerr << "\t-> Merging all reads into a single read group named: " << rgs[0] << " and sequencing cycles: " << cycles[0]  << '\n';
    settings.args_stream << "\t-> Merging all reads into a single read group named: " << rgs[0] << " and sequencing cycles: " << cycles[0]  << '\n';
  }

  // go to the correct chromosome
  // https://github.com/gatoravi/bam-parser-tutorial/blob/master/parse_bam.cc
  hts_itr_t *iter = sam_itr_querys(idx, header, settings.chrom.c_str());

  // load ref of that chromosome
  faidx_t *fai = ref_init(settings.reference_fn);

  char *ref = fetch_chrom(fai, settings.chrom);

  size_t seq_len = faidx_seq_len(fai, settings.chrom.c_str());

  if (!settings.exclude_bed_fn.empty()) {
    std::cerr << "\t-> Masking genomic BED regions (-e) " << settings.exclude_bed_fn << '\n';
    settings.args_stream << "\t-> Masking genomic BED regions (-e) " << settings.exclude_bed_fn << '\n';
    filter_ref_bed(settings.chrom, settings.exclude_bed_fn, ref);
  }

  if (!settings.exclude_sites_fn.empty()) {
    std::cerr << "\t-> Masking genomic sites sites (-E) " << settings.exclude_sites_fn << '\n';
    settings.args_stream << "\t-> Masking genomic sites sites (-E) " << settings.exclude_sites_fn << '\n';
    filter_ref_sites(settings.chrom, settings.exclude_sites_fn, ref);
  }

  std::vector<std::vector<double>> tmf;
  tmf.resize(rgs.size());
  std::vector<std::vector<int>> tm;
  tm.resize(rgs.size());
  for (size_t i=0; i<rgs.size(); i++){
    tm[i] = init_tallymat<int>(settings.max_pos_to_end);
  }

  my_cov_rg cov_rg(rgs.size());

  const keeplist_map cpg_map = get_cpg_chrom_pos(ref, seq_len);
  const std::vector<int> cpg_bool = get_cpg_chrom_bool(ref, seq_len);
  std::cerr << "\t-> " << cpg_map.size() << " CpG's in chrom: " << settings.chrom << '\n';
  settings.args_stream << "\t-> " << cpg_map.size() << " CpG's in chrom: " << settings.chrom << '\n';

  // do not included CnonCpGs if deamrates are already provided.
  std::vector<size_t> exclude_CnonCpGs(rgs.size(), 0);
  for (size_t i=0; i<rgs.size(); i++){
    if((!settings.deamrates_filename.empty()) && check_file_exists(settings.deamrates_filename)){
      exclude_CnonCpGs[i] = 1;
      std::cerr << "\t-> No need to include CnonCpGs for " << rgs[i] << " as -D is provided: " << settings.deamrates_filename << '\n';
      settings.args_stream << "\t-> No need to include CnonCpGs for " << rgs[i] << " as -D is provided: " << settings.deamrates_filename << '\n';

    } else if (check_file_exists(settings.outbase+"."+rgs[i]+".deamrates")){
      exclude_CnonCpGs[i] = 1;
      std::cerr << "\t-> No need to include CnonCpGs for " << rgs[i] << " as " << settings.outbase+"."+rgs[i]+".deamrates" << " exists" << '\n';
      settings.args_stream << "\t-> No need to include CnonCpGs for " << rgs[i] << " as " << settings.outbase+"."+rgs[i]+".deamrates" << " exists" << '\n';
    }
  }


  std::vector<per_site> data;
  data.resize(cpg_map.size());
  size_t nocpg_to_include = 1e6;
  per_site_nocpg nocpg_data(nocpg_to_include);

  bam1_t *rd = bam_init1();
  int reads = 1e6;
  int counter = 0;
  alignment_data d;
  size_t trashed=0, reads_skipped=0;
  size_t startpos, endpos;
  bool skipread=true;
  time_t start_time_load_data,end_time_load_data;

  uint8_t *rgptr;
  std::string rgname;
  size_t rgname_idx;
  std::string rname;
  time(&start_time_load_data);
  while (sam_itr_next(in, iter, rd) >= 0) {
    skipread = true;
    if (counter % reads == 0 && reads) {
      std::cerr << "\t-> " << counter << " reads processed and " << trashed << " discarded. " << '\r';
    }
    counter++;

    // this is 50% faster than analyzing every single read.
    startpos = rd->core.pos;
    endpos = bam_endpos(rd);

    // auto hit = std::find(cpg_bool.begin()+startpos, cpg_bool.begin()+endpos, 1);
    // if(hit!=cpg_bool.begin()+endpos){
    //   skipread=false;
    // }
    for (size_t i=startpos; i<=endpos; i++){
      if(cpg_bool[i]==1){
        skipread=false;
        break;
      }
    }

    if(skipread){
      reads_skipped++;
      continue;
    }

    if (rd->core.l_qseq < settings.minreadlength ||
        rd->core.qual < settings.minmapQ ||
        (rd->core.flag & settings.flags_off) != 0 ||
        rd->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) {
      trashed++;
      continue;
    }

    // when adding the rd to data, we should also add the RG as a vector. It
    // could be integers referecing the all the RG found in the header..

    if (rg_split) {
      rgptr = bam_aux_get(rd, "RG");
      if (rgptr == NULL) {
        // in the case norg exists. dont know if i should trash or put in
        // UNKWOWN
        rname = std::string(bam_get_qname(rd));
        std::cerr << "[No_found_RG] :: " << rname << '\n';
        trashed++;
        continue;
      } else {
        // moves past the Z by +1. all
        // colons are gone already
        rgname = std::string((const char *)(rgptr + 1));
        auto match = std::find(rgs.begin(), rgs.end(), rgname);
        if (match != rgs.end()) {
          rgname_idx = distance(rgs.begin(), match);
        } else { // in the case it cannot find the rg. The read will be trashed
          trashed++;
          continue;
        }
      }
   } else {
     rgname=ALL_RG;
     rgname_idx=ALL_DEAMMETH_RG_IDX;
   }

   d = align_read(rd, ref);
   add_aligned_data(settings, d, cpg_map, ref, data, tm[rgname_idx], nocpg_data, cov_rg, rgname_idx, cycles[rgname_idx], exclude_CnonCpGs[rgname_idx]);
  }
  time(&end_time_load_data);

  std::cerr << "\t-> Processed: " << counter  << ". Reads filtered: " <<
    trashed << ". Reads skipped (nocpg overlap): " << reads_skipped <<
    ". Loaded in " << difftime(end_time_load_data, start_time_load_data) << " seconds." << '\n';

  settings.args_stream << "\t-> Processed: " << counter  << ". Reads filtered: " <<
    trashed << ". Reads skipped (nocpg overlap): " << reads_skipped <<
    ". Loaded in " << difftime(end_time_load_data, start_time_load_data) << " seconds." << '\n';
  for (size_t i=0; i<rgs.size(); i++){
    std::cerr << "\t-> Total_Observations RG: "<< rgs[i] << " CpG: " << cov_rg.cpg[i] << " CpGCoverage: " << (double)cov_rg.cpg[i]/(double)cpg_map.size() <<  ". CnonCpG: " << cov_rg.nocpg[i] << '\n';
    settings.args_stream << "\t-> Total_Observations RG: "<< rgs[i] << " CpG: " << cov_rg.cpg[i] << " CpGCoverage: " << (double)cov_rg.cpg[i]/(double)cpg_map.size() <<  ". CnonCpG: " << cov_rg.nocpg[i] << '\n';
  }
  for (size_t i=0; i<rgs.size(); i++){
    std::cerr << "\t-> Dumping count file: " << settings.outbase << "." << rgs[i] <<  ".tallycounts" << '\n';
    settings.args_stream << "\t-> Dumping count file: " << settings.outbase << "." << rgs[i] <<  ".tallycounts" << '\n';
    dump_count_file(settings, tm[i], rgs[i]);
    tmf[i] = read_count_file(settings, rgs[i]);
  }

  // mark_no_deam_CT_GA(settings, data);

  std::vector<std::vector<double>> param_deam;
  param_deam.resize(rgs.size());

  settings.args_stream << std::flush;
  for (size_t i=0; i<rgs.size(); i++){
    if((!settings.deamrates_filename.empty()) && check_file_exists(settings.deamrates_filename)){
      std::cerr << "\t-> Loading deamination rates from " << settings.deamrates_filename << '\n';
      std::cerr << "\t-> Make sure that the file contains the same number of pos to include. Deammeth does not check that" << '\n';
      settings.args_stream << "\t-> Loading deamination rates from " << settings.deamrates_filename << '\n';
      param_deam[i] = load_deamrates_f(settings);
    } else if (check_file_exists(settings.outbase+"."+rgs[i]+".deamrates")){
      std::cerr << "\t-> Loading deamination rates from " << settings.outbase+"."+rgs[i]+".deamrates" << '\n';
      settings.args_stream << "\t-> Loading deamination rates from " << settings.outbase+"."+rgs[i]+".deamrates" << '\n';
      param_deam[i] = load_deamrates(settings, rgs[i]);
    }else {
      std::cerr << "\t-> Starting Optim of deamination rates. RG: " << rgs[i] << '\n';
      settings.args_stream << "\t-> Starting Optim of deamination rates. RG: " << rgs[i] << '\n';
      run_deamrates_optim(settings, tmf[i], data, nocpg_data, rgs[i], i, cov_rg);
      std::cerr << "\t-> Dumping deamination parameters to " << settings.outbase+"."+rgs[i]+".deamrates" << '\n';
      settings.args_stream << "\t-> Dumping deamination parameters to " << settings.outbase+"."+rgs[i]+".deamrates" << '\n';
      param_deam[i] = load_deamrates(settings, rgs[i]);
    }
  }
  std::vector<pre_calc_per_site> mle_data;

  double deamin_methylated, deamin_unmethylated;
  double noM, M;
  std::cerr << "\t-> Precalculating many things for ML of F." << '\n';
  mle_data.reserve(data.size());
  for (const auto &site : data) {
    // FIXME: do know if we should exclude or include missing sites.
    // For window sizes, it doesn't matter
    // for ncpg it means that we do not always analyze the same cpgs per sample.
    if(!site.depth){
      continue;
    }
    pre_calc_per_site d(site);
    double  prob_g1, prob_g2;
    auto g12 = SEVEN_DINUCL_GENOTYPES.begin();
    for (size_t i = 0; i < site.depth; i++) {
      deamin_unmethylated = param_deam[site.rgs[i]][get_param_idx(settings.max_pos_to_end, UNMETHSTATE, site.pos_to_end[i], site.prime[i])];
      // 1 should be replaced by do_haploid_model
      if(do_haploid_model){
        noM = calc_prob_obs_base(site.seqerrors[i], site.base_compos[i],
                                 deamin_unmethylated);
      } else {
        prob_g1 =
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].first,
                      site.bases[i].first, g12->first.first) +
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].first,
                      site.bases[i].first, g12->first.second);
        prob_g2 =
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].second,
                      site.bases[i].second, g12->second.first) +
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].second,
                      site.bases[i].second, g12->second.second);

        noM = prob_g1 * prob_g2;
      }
      d.pre_noM.push_back(noM);

      deamin_methylated = param_deam[site.rgs[i]][get_param_idx(settings.max_pos_to_end, METHSTATE, site.pos_to_end[i], site.prime[i])];

      // 1 should be replaced by do_haploid_model
      if(do_haploid_model){
        M = calc_prob_obs_base(site.seqerrors[i], site.base_compos[i], deamin_methylated);
      } else {
        prob_g1 =
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_methylated, site.quals[i].first,
                      site.bases[i].first, g12->first.first) +
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_methylated, site.quals[i].first,
                      site.bases[i].first, g12->first.second);
        prob_g2 =
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_methylated, site.quals[i].second,
                      site.bases[i].second, g12->second.first) +
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_methylated, site.quals[i].second,
                      site.bases[i].second, g12->second.second);

        M = prob_g1 * prob_g2;
      }
      d.pre_M.push_back(M);

      d.maperrors.push_back(site.maperrors[i]);
    }
    for (auto it=SEVEN_DINUCL_GENOTYPES.begin()+1; it!=SEVEN_DINUCL_GENOTYPES.end(); it++){
      double res=0;
      for (size_t i = 0; i < site.depth; i++) {
        deamin_unmethylated =
            param_deam[site.rgs[i]]
                      [get_param_idx(settings.max_pos_to_end, UNMETHSTATE,
                                     site.pos_to_end[i], site.prime[i])];
        prob_g1 =
            // std::log((1 - site.maperrors[i]) *
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].first,
                      site.bases[i].first, it->first.first) +
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].first,
                      site.bases[i].first, it->first.second);
        // + site.maperrors[i] * DINUCL_FLAT_PRIOR; //);
        prob_g2 =
            // std::log((1 - site.maperrors[i]) *
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].second,
                      site.bases[i].second, it->second.first) +
            0.5 * base_condition_one_genotype(
                      site.strand[i], deamin_unmethylated, site.quals[i].second,
                      site.bases[i].second,
                      it->second.second);
        // + site.maperrors[i] * DINUCL_FLAT_PRIOR);
        res += std::log((1-site.maperrors[i]) * (prob_g1 * prob_g2) +
                        site.maperrors[i] * DINUCL_FLAT_PRIOR);
      }
      size_t idx_prior = std::distance(SEVEN_DINUCL_GENOTYPES.begin(), it);
      d.remaining_dinucl_genotypes.push_back(res + LOG_PRIORS[idx_prior]);
    }
    mle_data.push_back(d);
    // if(d.position==112256){
    //   print_d(d);
    // }
  }
  data.clear();
  if(settings.bed_f.empty()){
    std::cerr << "\t-> Dumping MLE of F to " << settings.outbase << ".F" << '\n';
    settings.args_stream << "\t-> Dumping MLE of F to " << settings.outbase << ".F" << '\n';
    run_mle(settings, mle_data);

  } else {
    std::cerr << "\t-> Dumping MLE of F to " << settings.outbase << ".BED.F" << '\n';
    settings.args_stream << "\t-> Dumping MLE of F to " << settings.outbase << ".BED.F" << '\n';
    run_mle_bed(settings, mle_data, bed_coord);
  }
  std::cerr << "\t-> Cleaning up." << '\n';
  settings.args_stream << std::flush;
  mle_data.clear();
  free(ref);
  hts_itr_destroy(iter);
  hts_idx_destroy(idx);
  bam_destroy1(rd);
  bam_hdr_destroy(header);
  sam_close(in);
  return 1;
}


void mark_no_deam_CT_GA(general_settings & settings, std::vector<per_site> & data){
  std::default_random_engine generator;
  std::uniform_real_distribution<double> distribution(0.0,1.0);

  std::string filename = settings.outbase +".excludedCpGs";
  std::ofstream f (filename.c_str());
  checkfilehandle(f, filename);


  for (auto & site: data){
    double prob_seqerrors=1;
    for (size_t i=0; i<site.depth; i++){
      // + strand and g->a
      if(site.strand[i]==0 && site.bases[i].second==0){
        site.no_deam_CT_GA++;
        prob_seqerrors *= site.seqerrors[i]/3.0;
      }
      // - strand and c->t
      if(site.strand[i]==1 && site.bases[i].first==3){
        site.no_deam_CT_GA++;
        prob_seqerrors *= site.seqerrors[i]/3.0;
      }
    }
    // we keep a site with probability prop_seqerrors i.e the observations are a series of sequencing errors and not true variants.
    double accept_val = distribution(generator);
    if(site.no_deam_CT_GA>1 && accept_val>prob_seqerrors){
      site.exclude_for_deam = true;
      f << settings.chrom << " " << site.position << '\n';
    }
  }
}










double objective_func_F_haploid(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
  F_void *d = static_cast<F_void*>(my_func_data);
  d->iteration++;
  double res, ll=0;
  double noM, M;
  double mape;
  if(!grad.empty()){
    grad[0] = 0;
  }
  pre_calc_per_site * site;
  for (const auto idx : d->mle_data->idx_to_include){
    site = &d->data->at(idx);
    for (size_t i=0; i < site->depth; i++){
      noM =  site->pre_noM[i];
      M = site->pre_M[i];
      mape = site->maperrors[i];
      res = (1-mape) * ((1-x[0]) * noM + x[0] * M) + mape*BASE_FREQ_FLAT_PRIOR;
      ll += std::log(res);
      // wolfram: derivative  ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // slope += ((1-mape) * (M - noM)) / ((mape-1) * ((f-1) * noM - f * M) + mape*BASE_FREQ_FLAT_PRIOR);
      // beloved wolfram got it right, just modified a few things for readability
      if(!grad.empty()){
        grad[0] += ((1-mape) * (M - noM)) / res;
      }
      // https://www.symbolab.com/
      // wolfram: second derivative ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // second derivative:
      // second_der += (std::pow((1-mape),2) * std::pow((M - noM),2)) / std::pow((-M * (mape-1) * f + noM * (mape-1) * (1-f) + mape*BASE_FREQ_FLAT_PRIOR), 2);
    }
  }
  return ll;
}

double objective_func_F_second_deriv_haploid(const double & f, std::vector<pre_calc_per_site> & data, per_mle_run &mle_data){
  double slope=0;
  double nominator, denominator;
  double noM, M;
  double mape;
  pre_calc_per_site * site;
  for (const auto idx : mle_data.idx_to_include){
    site = &data.at(idx);
    for (size_t i=0; i < site->depth; i++){
      noM =  site->pre_noM[i];
      M = site->pre_M[i];
      mape = site->maperrors[i];
      nominator = std::pow((-mape + 1 ), 2) * std::pow((M - noM ), 2);
      denominator = std::pow((1-mape) * (noM * (1-f) + f * M) + mape*BASE_FREQ_FLAT_PRIOR, 2);
      slope += -nominator/denominator;
    }
  }
  return slope;
}
