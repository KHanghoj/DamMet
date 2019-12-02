double loglikelihood_after_mle(const std::vector<double> &x, pre_calc_per_site * site){
  double log_like_cgcg_geno=0;
  double read_like ;
  double noM, M, mape ;
  for (size_t i=0; i < site->depth; i++){
    noM =  site->pre_noM[i];
    M = site->pre_M[i];
    mape = site->maperrors[i];
    read_like = (1-mape) * ((1-x[0]) * noM + x[0] * M) + mape*DINUCL_FLAT_PRIOR;
    log_like_cgcg_geno += std::log(read_like);
  }
  return log_like_cgcg_geno + LOG_PRIORS[0];
}

double objective_func_F_second_deriv(const double & f, std::vector<pre_calc_per_site> & data, per_mle_run &mle_data){
  double total_d2=0;
  // double nominator, denominator;
  double log_like_cgcg_geno;
  double reads_grad, geno_grad;
  double read_der, read_like, log_like_genos ;
  double noM, M;
  double mape;
  double numerator_read_der2, denominator_read_der2, reads_der2;
  double geno_grad_der2, like_genos, like_genos_der2;
  pre_calc_per_site * site;
  for (const auto & idx : mle_data.idx_to_include){
    site = &data.at(idx);
    log_like_cgcg_geno = 0;
    reads_grad=0;
    reads_der2=0;
    for (size_t i=0; i < site->depth; i++){
      noM =  site->pre_noM[i];
      M = site->pre_M[i];
      mape = site->maperrors[i];

      read_der = (1-mape) * (M - noM);
      read_like = (1-mape) * ((1-f) * noM + f * M) + mape*DINUCL_FLAT_PRIOR;
      reads_grad += read_der/read_like;

      numerator_read_der2 = std::pow((1-mape),2) * std::pow((M - noM),2);
      denominator_read_der2 = std::pow((1-mape) * (noM * (1-f) + f * M) + mape*DINUCL_FLAT_PRIOR, 2);
      reads_der2 += - (numerator_read_der2 / denominator_read_der2);
      // \frac{\partial }{\partial \:x^2}\left(ln\left(\left(1-w\right)\cdot \left(\left(1-x\right)\:\cdot \:n\:+\:\left(x\:\cdot \:M\right)\right)\:+\:w\cdot \:p\right)\:\:\right)

      log_like_cgcg_geno += std::log(read_like);
      // https://www.symbolab.com/
      // wolfram: second derivative ln((1-w)*((1-x) * k + (x * h)) + w*p)
      // second derivative:
      // second_der += (std::pow((1-mape),2) * std::pow((M - noM),2)) / std::pow((-M * (mape-1) * f + noM * (mape-1) * (1-f) + mape*DINUCL_FLAT_PRIOR), 2);
    }
    // exp(log_like_cgcg_geno + LOG_PRIORS) * d/df all_reads. this is only including cgcg as the derivatives of the remaining are 0 as they are unrelated to F.

    //g(x) first derivative     z(x)                              y(x)
    geno_grad = std::exp(log_like_cgcg_geno + LOG_PRIORS[0]) * reads_grad;

    // g'(x) (of the first derivative) -> z(x)*y'(x) + z'(x) * y(x)
    geno_grad_der2 = std::exp(log_like_cgcg_geno + LOG_PRIORS[0]) * reads_der2 + (geno_grad * reads_grad);
    log_like_genos = log_like_cgcg_geno + LOG_PRIORS[0];
    // summation of the likelihood of the remaining genotypes
    for (auto it=site->remaining_dinucl_genotypes.begin(); it!=site->remaining_dinucl_genotypes.end(); it++){
       log_like_genos = oplusnatl(log_like_genos, *it);
    }

    // f(x) first derivative is 1/below
    like_genos = std::exp(log_like_genos);
    // f'(x) of the first derivative
    like_genos_der2 = (- 1.0/(std::pow(like_genos,2))) * geno_grad;

    //             f'(x)           g(x)           f(x)              g'(x)
    total_d2 += like_genos_der2 * geno_grad + (1.0/like_genos) * geno_grad_der2;
  }
  return total_d2;
}


// void run_mle_bed(general_settings & settings,
//                  std::vector<pre_calc_per_site> & pre_calc_data,
//                  std::vector<std::pair<size_t, size_t>> & bed_coord) {
//   std::string filename = settings.outbase + ".BED.F";
//   std::ofstream f (filename.c_str());
//   checkfilehandle<std::ofstream>(f, filename);
//   // nlopt::opt opt(nlopt::LN_BOBYQA, 1);
//   nlopt::opt opt(nlopt::LD_MMA, 1);
//   std::vector<double> lb(1, SMALLTOLERANCE), up(1, 1-SMALLTOLERANCE);
//   opt.set_maxeval(1000);
//   opt.set_lower_bounds(lb);
//   opt.set_upper_bounds(up);
//   opt.set_xtol_abs(0);
//   opt.set_ftol_abs(1e-15);
//   opt.set_xtol_rel(0);
//   opt.set_ftol_rel(0);
//   for (const auto & bed : bed_coord){
//     std::vector<size_t> sites_to_include;
//     for (size_t curr_idx = 0; curr_idx < pre_calc_data.size(); curr_idx++) {

//       if (pre_calc_data[curr_idx].position >= bed.first && pre_calc_data[curr_idx].position < bed.second){

//         sites_to_include.push_back(curr_idx);
//       }
//       if(pre_calc_data[curr_idx].position >= bed.second){
//         break;
//       }

//     }
//     if( sites_to_include.size()==0){
//       f << settings.chrom << ":" << bed.first<<"-"<<bed.second << " " << "NO_SITES_IN_THE_REGION NOT_USED" << '\n';
//       continue;
//     }

//     // it is stupid to run over data again, but too tired to make new structs for the BED setup.
//     // FIXME: make a struct that does not need to be initialized, cause then it can be used in the for loop above.
//     per_mle_run mle_data(pre_calc_data[sites_to_include[0]], sites_to_include[0]);
//     for (auto it=sites_to_include.begin()+1; it!=sites_to_include.end();it++){
//       update_mle_data(mle_data, pre_calc_data[*it], *it);
//     }

//     if(mle_data.total_depth == 0){
//       f << settings.chrom << ":" << bed.first<<"-"<<bed.second << " " << "NO_DATA_IN_THE_REGION NOT_USED"  << '\n';
//       continue;
//     }

//     double minf;
//     F_void void_stuff(&settings, &pre_calc_data, &mle_data);
//     // std::vector<double> param (1, 0.5);
//     std::vector<double> param (1, SMALLTOLERANCE);
//     nlopt::result result;
//     double second_der, error;
//     size_t iterations;

//     // // temp
//     // std::vector<double> dummy;
//     // for (double i=0; i<=1; i+=0.01){
//     //   std::vector<double> x(i);
//     //   std::cout << i << " " << std::setprecision(10) << objective_func_F(param, dummy, &void_stuff) << std::endl;
//     // }
//     // loglikelihood_after_mle(param[0],
//     // // temp end

//     if(do_haploid_model){
//       opt.set_max_objective(objective_func_F_haploid, &void_stuff);
//       result = opt.optimize(param, minf);
//       second_der = objective_func_F_second_deriv_haploid(param[0], pre_calc_data, mle_data);
//     } else {
//       opt.set_max_objective(objective_func_F, &void_stuff);
//       result = opt.optimize(param, minf);
//       second_der = objective_func_F_second_deriv(param[0], pre_calc_data, mle_data);
//     }
//     error = 1.96/std::sqrt(-second_der);
//     iterations=void_stuff.iteration;

//     f << settings.chrom << ":" << bed.first<<"-"<<bed.second
//       << " N_CpGs: " << mle_data.n_cpgs
//       << " Depth: " << mle_data.total_depth
//       << " Distance: " << mle_data.max_pos - mle_data.min_pos
//       << " ll: " << minf
//       << " f: " << param[0]
//       << " f(95%conf): " << param[0]-error  << "," <<  param[0]+error
//       << " iterations: " << iterations
//       << " optim_return_code: " << result
//       << " Incl_pos: " << mle_data.positions[0];
//       for (auto i=mle_data.positions.begin()+1; i!=mle_data.positions.end(); i++){
//         f << ","<< *i;
//       }
//       f << std::endl;

//   }
//   f.close();
// }

double calc_prob_obs_base(const double & seqerror, const size_t & base_compos, const double & deamin_rate){
  double res;
  if(base_compos == 0){  // CC :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
    res = ((1-deamin_rate) * (1-seqerror)) + (deamin_rate * seqerror/3) + (seqerror * seqerror/3);
  }else if(base_compos == 1){  // CT :   seqerror/3 comes from AGT -> C. seqerror -> C -> ACG
    res = (deamin_rate * (1-seqerror)) + ((1-deamin_rate) * seqerror/3) + (seqerror * seqerror/3);
  } else { // mutation
    res = ((1-deamin_rate) * seqerror/3) + (deamin_rate * seqerror/3) + (seqerror * seqerror/3);
  }
  return res;
}



// const std::vector<int> get_cpg_chrom_pos (const char * ref, const size_t & seq_len){
//   size_t counter = 0;
//   for (size_t i=0; i<seq_len-1; i++){
//     if (refToInt[(int)ref[i]] == 1 && refToInt[(int)ref[i+1]] == 2){
//       res.push_back(i);
//     }
//   }
//   return res;
// }

// const v_un_ch get_c_and_cpg_chrom(const char * ref, const size_t & seq_len){
//   v_un_ch res(seq_len, 0);
//   for (size_t i=0; i<seq_len-1; i++){
//     if (refToInt[(int)ref[i]] == 1){

//       if (refToInt[(int)ref[i+1]] == 0){
//         res[i] = 1;
//         res[i+1] = 1;
//       } else if (refToInt[(int)ref[i+1]] == 1){
//         res[i] = 2;
//         res[i+1] = 2;
//       } else if (refToInt[(int)ref[i+1]] == 2){
//         res[i] = 3;
//         res[i+1] = 3
//       } else if (refToInt[(int)ref[i+1]] == 3){

//       }
//     }
//   }
//   return res;
// }



// std::vector<int> idx_to_params(const int & idx){
//   std::div_t p = std::div(idx, READPOS_MULT);
//   std::div_t pr = std::div(p.rem, PRIME_MULT);
//   std::div_t strand = std::div(pr.rem, STRAND_MULT);
//   std::div_t B1 = std::div(strand.rem, B1_MULT);
//   std::div_t B2 = std::div(B1.rem, B2_MULT);
//   std::div_t B3 = std::div(B2.rem, B3_MULT);
//   return std::vector<int> {p.quot, pr.quot, strand.quot, B1.quot, B2.quot, B3.quot, B3.rem};
// }
