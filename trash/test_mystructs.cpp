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

void r1(const size_t & nsites, const size_t & nobs ){
  std::vector<std::unique_ptr<Site>> all;
  for (size_t pos=0; pos<nsites; pos++){
    all.push_back(std::make_unique<Site>(pos, 1));
    for (size_t i=0; i<nobs; i++){
      // all.back()->data.push_back(std::make_unique<Obs>(0, 1, 10, 1, 1, 40, 30, 1));
      all.back()->data.push_back(std::make_unique<Obs>(0, 1, 0, 0, 0, 0, 0, 0));
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
  std::cerr << "sleeping" << '\n';
  std::this_thread::sleep_for(std::chrono::milliseconds(10 * s));
}

void r2(const size_t & nsites, const size_t & nobs ){
  std::vector<std::unique_ptr<Site>> all;
  Site d(0,1);
  for (size_t pos=0; pos<nsites*nobs; pos++){
    d.data.push_back(std::make_unique<Obs>(0, 1, 0, 0, 0, 0, 0, 0));
  }
  double mem;
  int s=1000;
  std::cerr << d.data.size() << '\n';
  mem = sizeof(std::unique_ptr<Obs>) * d.data.size();
  mem /= MEGABYTE;
  std::cerr << mem << '\n';
  std::cerr << "sleeping" << '\n';
  std::this_thread::sleep_for(std::chrono::milliseconds(10 * s));
}

int main(int argc, char * argv[]){
  if (argc != 3){
    std::cerr << "need 'nsites' and 'nobs' argument only" << '\n';
    exit(EXIT_FAILURE);
}
  size_t nsites = atoi(argv[1]);
  size_t nobs = atoi(argv[2]);
  r1(nsites, nobs);
  r2(nsites, nobs);

}
