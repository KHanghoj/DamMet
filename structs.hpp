// NEW structs

using unint = unsigned int;

/// STRUCTS
struct my_cov_rg {
  size_t nocpg=0, cpg=0;
} ;

struct alignment_data {
  alignment_data () {
    t_seq.reserve(100), t_posi.reserve(100), t_isop.reserve(100), t_qs.reserve(100), t_ref.reserve(100), t_positions.reserve(100);
  }
  unint strand, mapQ, n_nucleotides;
  std::vector<unint> t_seq, t_posi, t_isop, t_qs, t_ref, t_positions;

};



struct ObsF {
  ObsF(const double &_M,
       const double &_noM,
       const double &_me
       ):
    M(_M), noM(_noM), me(_me){}
  double M, noM;
  double me;
};

struct Obs {
  Obs(const unint &_pr, // prime 0,1
      const unint &_st, // strand 0,1
      const unint &_rp, // readpos int
      const unint &_rl,  // read length 0,1
      const unint &_bc, // base composition 0,1,2
      const unint &_se, // seq error 1-40
      const unint &_me // mapping error 1-40
      ):
    pr(_pr), st(_st),
    rp(_rp), rl(_rl),
    bc(_bc), se(_se),
    me(_me) {

  }

  unint pr : 1;  // 0..1  (1 bits)
  unint st : 1;  // 0..1  (1 bits)
  unint rp : 6;  // 0..64 (6 bits) 8
  unint rl : 1;  // 0..1  (1 bits)
  unint bc : 2;  // 0,1,2 (2 bits) 3
  unint se : 6;  // 0-40  (6 bits)
  unint me : 6;  // 0-40  (6 bits)
};

using uni_ptr_obs = std::vector<std::unique_ptr<Obs>>;

struct Site_s {
  Site_s(){
    remaining_dinucl_genotypes.resize(6, 0);
  };

  void load_data(const double &_M, const double &_noM, const double &_me){
    data.push_back(ObsF(_M, _noM, _me));
    depth++;
  };
  size_t pos, depth=0;
  std::vector<double> remaining_dinucl_genotypes;
  std::vector<ObsF> data;

};

struct per_mle_run {
  std::vector<std::unique_ptr<Site_s>> to_include; // , positions;
  std::vector<size_t> positions;
  Site_s center_site;
  size_t total_depth, min_pos, max_pos, curr_pos, n_cpgs;
  per_mle_run (const Site_s & _d){
    n_cpgs = 1;
    to_include.emplace_back(std::make_unique<Site_s>(_d));
    center_site = _d;
    positions.push_back(_d.pos);
    total_depth = _d.depth;
    min_pos = _d.pos;
    max_pos = _d.pos;
    curr_pos = _d.pos;
  };

  void update_left(const Site_s & _d){
    min_pos = _d.pos;
    stats_update(_d);
  }

  void update_right(const Site_s & _d){
    max_pos = _d.pos;
    stats_update(_d);
  }

  void stats_update(const Site_s & _d){
    to_include.emplace_back(std::make_unique<Site_s>(_d));
    positions.push_back(_d.pos);
    total_depth += _d.depth;
    n_cpgs++;
  }

} ;

struct rgs_info {
  bool rg_split=false;
  std::vector<std::string> rgs;
  std::vector<size_t> cycles;
  size_t n = 0;
};

struct deamrates_void {
  deamrates_void(general_settings * _settings,
                 uni_ptr_obs &_cpg_data,
                 uni_ptr_obs &_nocpg_data){
    settings = _settings;
    cpg_data = std::move(_cpg_data);
    nocpg_data = std::move(_nocpg_data);
  };
  general_settings * settings;
  uni_ptr_obs cpg_data, nocpg_data;
  size_t iteration=0;

};


struct F_void {
  general_settings * settings;
  std::vector<std::unique_ptr<Site_s>> to_include;
  size_t iteration;
  F_void(general_settings * _s, std::vector<std::unique_ptr<Site_s>> & _m){
    settings = _s;
    to_include = std::move(_m);
    iteration = 0;
  }

  F_void(std::vector<std::unique_ptr<Site_s>> & _m){
    to_include = std::move(_m);
  }
};


// multi thread
struct job_deamrates {

  job_deamrates(rgs_info &_rgs){
    rgs = _rgs;
    tm.resize(rgs.n);
    cpg_data.resize(rgs.n);
    nocpg_data.resize(rgs.n);
    cov_rg.resize(rgs.n);
  }

  void add_chrom(std::string &c){
    chroms.push_back(c);
  };
  rgs_info rgs;
  std::vector<std::vector<int>> tm;
  std::vector<uni_ptr_obs> cpg_data, nocpg_data;
  std::vector<my_cov_rg> cov_rg;
  std::vector<std::string> chroms;
};
