#include <iostream> // stdout/stdin/stderr
#include <memory>
#include <vector>
// bit_fields1.cpp
// compile with: /LD
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
    data.reserve(20);
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

int main(){
  Site d(1000, 1);
  for (size_t i=0; i<20; i++){
    d.data.push_back(std::make_unique<Obs>(0, 1, 10, 1, 1, 40, 30, 1));
    d.depth++;
  }
  std::cerr << d.depth++ << '\n';
  for (auto &x: d.data){
    std::cerr << d.pos << " " << d.ref << " ";
    print(x);

  }

}
