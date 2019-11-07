#include <iostream> // stdout/stdin/stderr
#include <memory>
#include <vector>
#include <chrono>
#include <thread>
// bit_fields1.cpp
// compile with: /LD
const unsigned long MEGABYTE = 1024 * 1024;


using unint = unsigned int;

struct Obs {
  Obs(const unint &_pr, // prime 0,1
      const unint &_st, // strand 0,1
      const unint &_rp, // readpos int
      const unint &_rl,  // read length 0,1
      const unint &_bc, // base composition 0,1,2
      const unint &_se, // seq error 1-40
      const unint &_me, // mapping error 1-40
      const unint &_rg  // read group idx 1-x
      ):
    pr(_pr), st(_st),
    rp(_rp), rl(_rl),
    bc(_bc), se(_se),
    me(_me), rg(_rg) {}

  unint pr : 1;  // 0..1  (1 bits)
  unint st : 1;  // 0..1  (1 bits)
  unint rp : 6;  // 0..64 (6 bits) 8
  unint rl : 1;  // 0..1  (1 bits)
  unint bc : 2;  // 0,1,2 (2 bits) 3
  unint se : 6;  // 0-40  (6 bits)
  unint me : 6;  // 0-40  (6 bits)
  unint rg : 6;  // 0-40  (6 bits) 18
};

struct Site {
  Site(const unint &_pos,
       const unint &_ref):
    pos(_pos), ref(_ref){
    // data.reserve(20);
  }
  int pos, depth=0;
  unint ref : 3;  // 0-5 A,C,G,T,N,Del
  std::vector<std::unique_ptr<Obs>> data;
};

struct Siteold {
  Siteold(const unint &_pos,
          const unint &_ref):
    pos(_pos), ref(_ref){  }
  int pos, depth=0, ref ;
  std::vector<int> pr, st, rp, rl, bc, se, me, rg;
};


void Siteold_update(Siteold & d,
                    const int &_pr, // prime 0,1
                    const int &_st, // strand 0,1
                    const int &_rp, // readpos int
                    const int &_rl,  // read length 0,1
                    const int &_bc, // base composition 0,1,2
                    const int &_se, // seq error 1-40
                    const int &_me, // mapping error 1-40
                    const int &_rg){  // read group idx 1-x):
  d.pr.push_back(_pr);
    d.st.push_back(_st);
    d.rp.push_back(_rp);
    d.rl.push_back(_rl);
    d.bc.push_back(_bc);
    d.se.push_back(_se);
    d.me.push_back(_me);
    d.rg.push_back(_rg);
    // data.reserve(20);
}

void print(std::unique_ptr<Obs> & d){
  std::cerr << d->pr << " ";
  std::cerr << d->se << " ";
  std::cerr << d->rp << " ";
  std::cerr << d->rl << " ";
  std::cerr << d->bc << " ";
  std::cerr << d->se << " ";
  std::cerr << d->me << " ";
  std::cerr << d->rg << '\n';
}

int r1(const size_t & nsites, const size_t & nobs ){
  std::vector<std::unique_ptr<Site>> all;
  for (size_t pos=0; pos<nsites; pos++){
    all.emplace_back(std::make_unique<Site>(pos, 1));
    for (size_t i=0; i<nobs; i++){
      all.back()->data.emplace_back(std::make_unique<Obs>(0, 1, 0, 0, 0, 0, 0, 0));
      all.back()->depth++;
    }
  }
  int s = 1000;
  double mem;
  std::cerr << all.size() << " " << all[0]->data.size() << '\n';
  mem = sizeof(std::unique_ptr<Site>) * all.size();
  mem *= sizeof(std::unique_ptr<Obs>) * all[0]->data.size();
  mem /= MEGABYTE;

  std::cerr << all[0]->depth++ << " " << mem << '\n';
  // std::cerr << "sleeping" << '\n';
  // std::this_thread::sleep_for(std::chrono::milliseconds(10 * s));
  return 1;
}

int r2(const size_t & nsites, const size_t & nobs ){

  // std::unique_ptr<Site> d(std::make_unique<Site>(0,1));
  std::vector<std::unique_ptr<Obs>> d;
  for (size_t pos=0; pos<nsites*nobs; pos++){
    d.emplace_back(std::make_unique<Obs>(0, 1, 10, 1, 1, 40, 30, 1));
  }
  double mem;
  int s=1000;
  std::cerr << d.size() << '\n';
  mem = sizeof(std::unique_ptr<Obs>) * (d.size());
  mem /= MEGABYTE;
  std::cerr << mem << '\n';
  // std::cerr << "sleeping" << '\n';
  // std::this_thread::sleep_for(std::chrono::milliseconds(10 * s));
  return 1;
}

int old(const size_t & nsites, const size_t & nobs){
  Siteold d(1,0);
  for (size_t pos=0; pos<nsites*nobs; pos++){
    Siteold_update(d, 0, 1, 10, 1, 1, 40, 30, 1);
  }
  return 1;
}

int main(int argc, char * argv[]){
  if (argc != 3){
    std::cerr << "need 'nsites' and 'nobs' argument only" << '\n';
    exit(EXIT_FAILURE);
}
  size_t nsites = atoi(argv[1]);
  size_t nobs = atoi(argv[2]);
  int a = r1(nsites, nobs);
  //std::cerr << "releasing memory" << '\n';
  //std::this_thread::sleep_for(std::chrono::milliseconds(5 * 1000));
  int b = r2(nsites, nobs);
  //std::cerr << "releasing memory2" << '\n';
  //std::this_thread::sleep_for(std::chrono::milliseconds(5 * 1000));
  int c = old(nsites, nobs);
  //std::cerr << "releasing memory2" << '\n';
  //std::this_thread::sleep_for(std::chrono::milliseconds(5 * 1000));
}



int main2(int argc, char * argv[]){
  if (argc != 3){
    std::cerr << "need 'nsites' and 'nobs' argument only" << '\n';
    exit(EXIT_FAILURE);
}
  size_t nsites = atoi(argv[1]);
  size_t nobs = atoi(argv[2]);

  std::vector<std::unique_ptr<Obs>> d;
  for (size_t pos=0; pos<nsites*nobs; pos++){
    d.emplace_back(std::make_unique<Obs>(0, 1, 10, 1, 1, 40, 30, 1));
  }
  std::cerr << "not releasing memory2" << '\n';
  std::this_thread::sleep_for(std::chrono::milliseconds(5 * 1000));
  std::vector<std::unique_ptr<Site>> all;
  for (size_t pos=0; pos<nsites; pos++){
    all.emplace_back(std::make_unique<Site>(pos, 1));
    for (size_t i=0; i<nobs; i++){
      all.back()->data.emplace_back(std::make_unique<Obs>(0, 1, 0, 0, 0, 0, 0, 0));
      all.back()->depth++;
    }
  }
  std::cerr << "not releasing memory" << '\n';
  std::this_thread::sleep_for(std::chrono::milliseconds(5 * 1000));
}
