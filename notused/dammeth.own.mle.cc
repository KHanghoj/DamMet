void run_mle(const general_settings & settings,
             const std::vector<pre_calc_per_site> & pre_calc_data) {

  for (size_t curr_idx = 0; curr_idx < pre_calc_data.size(); curr_idx++) {

    per_mle_run mle_data(pre_calc_data[curr_idx], curr_idx);
    size_t pos_to_left = (curr_idx > 0) ? curr_idx - 1 : 0;
    size_t pos_to_right = curr_idx + 1;
    bool can_go_left = false, can_go_right = false;
    // keep adding site to mle_data until we fulfill the requirements
    // while ((1 - mle_data.prob_no_deamin_event_meth) < settings.min_prob_deamin_event) {
    // MAX CPG's
    while (mle_data.n_cpgs < 20) {
      can_go_right = false;
      can_go_left = false;

      if (pos_to_right < pre_calc_data.size() - 1) {
        can_go_right = true;
      }

      if (pos_to_left >= 1) {
        can_go_left = true;
      }

      if (!pos_to_right && !pos_to_left) {
        // break ;; cannot extend anymore
        break;

      } else if (pos_to_right && !pos_to_left) {
        // take the one to the right
        update_mle_data(mle_data, pre_calc_data[pos_to_right], pos_to_right);
        pos_to_right++;
      } else if (!pos_to_right && pos_to_left) {
        // take the one to the left
        update_mle_data(mle_data, pre_calc_data[pos_to_left], pos_to_left);
        pos_to_left--;
      } else if (pos_to_right && pos_to_left) {
        // take the closest
        if (pre_calc_data[pos_to_right].position - mle_data.curr_pos <=
            mle_data.curr_pos - pre_calc_data[pos_to_left].position) {
          // take the one to the right
          update_mle_data(mle_data, pre_calc_data[pos_to_right], pos_to_right);
          pos_to_right++;
        } else {
          // take the one to the left
          update_mle_data(mle_data, pre_calc_data[pos_to_left], pos_to_left);
          pos_to_left--;
        }
      } else {
        // we are screwed. should not be possible
      }

    } // end while loop
    double ALPHA = 1e-8;
    double tole = 1e-10;
    double BETA = 0.5;

    // here starts the mle estimator.
    double f=0.5, fnew=0.5; // initial guess
    // fprintf(stderr, "startinglikeli: %f\n", oldll);
    size_t niter = 20000000;
    size_t c_iter=0;
    mle_ll_run ll_run, ll_runold;
    mle_loglike(ll_runold, f, pre_calc_data, mle_data);
    fnew = f + ALPHA * ll_runold.slope;
    fprintf(stderr, "startlik: %f f: %f fnew:%f\n", ll_runold.ll, f, fnew);
    double z = ll_runold.slope;
    for (f = 0; f<=1; f+=0.01){
    // for (c_iter = 0 ; c_iter < niter; c_iter++) {
    //   mle_loglike(ll_run, fnew, pre_calc_data, mle_data);
    //   z = BETA * z + ll_run.slope;
    //   fnew = f + ALPHA * z;
      BETA=0;
      mle_loglike(ll_run, f, pre_calc_data, mle_data);
      z = BETA * z + ll_run.slope;
      fnew = f + ALPHA * z;
      fprintf(stdout, "lik=%f diff=%e f=%f fnew=%f slope:%f second_der:%f\n", ll_run.ll, fabs(ll_run.ll - ll_runold.ll), f, fnew, ll_run.slope, ll_run.second_der);
      ll_runold = ll_run;

      // if(c_iter%10000==0 && c_iter)
      //   fprintf(stderr, "[%lu] lik=%f diff=%e f=%f slope:%f\n", c_iter, ll_run.ll, fabs(ll_run.ll - ll_runold.ll),
      //           fnew, ll_run.slope);
      // if (fabs(ll_run.ll - ll_runold.ll) < tole) {
      //   // fprintf(stderr,"breaking\n");
      //   break;
      // }
      // if (fnew>1){
      //   fnew = 1-tole;
      // }
      // if(fnew<0){
      //   fnew = tole;
      // }

      // if (ll_run.ll > ll_runold.ll){
      //   ALPHA *= .9;
      // }

      // f = fnew;
      // ll_runold = ll_run;
    }
    mle_ll_run ll_low, ll_high;
    double low = 0.0000001;
    double high = 0.9999999;
    mle_loglike(ll_low, low, pre_calc_data, mle_data);
    mle_loglike(ll_high, high, pre_calc_data, mle_data);
    std::cout << "Starting EM (contig: " << settings.chrom
              << ") with center position: " << mle_data.curr_pos
              << " ::: N_CpGs: " << mle_data.n_cpgs
              << " ::: Depth: " << mle_data.total_depth
              << " (CG:" << mle_data.summary[0]
              << ",TG:" << mle_data.summary[1]
              << ",{A,G}G:" << mle_data.summary[2]
              << ") ::: Prob_of_deam_meth: " << 1 - mle_data.prob_no_deamin_event_meth
              << " ::: Prob_of_deam_unmeth: " << 1 - mle_data.prob_no_deamin_event_unmeth
              << " ::: Min position: " << mle_data.min_pos
              << " ::: Max position: " << mle_data.max_pos
              << " ::: Distance (bp): " << mle_data.max_pos - mle_data.min_pos
              << " ::: ll: " << ll_run.ll
              << " ::: f: " << f
              << " ::: iterations: " << c_iter
              << " ::: ll0: " << ll_low.ll
              << " ::: ll1: " << ll_high.ll
              << '\n';
    std::cerr << "BREAKING" << '\n';
    break;
  }
}
