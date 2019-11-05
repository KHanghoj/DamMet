#include <iostream> // stdout/stdin/stderr
#include <memory>
#include <vector>
// bit_fields1.cpp
// compile with: /LD
using unint = unsigned int;

struct Date {
  Date(const unint &pnWD, const unint &pnMD,
       const unint &pnM, const unint &pnY):
    nWD(pnWD), nMD(pnMD), nM(pnM), nY(pnY) {}

  unint nWD : 3;    // 0..7   (3 bits)
  unint nMD : 6;    // 0..31  (6 bits)
  unint nM  : 5;    // 0..12  (5 bits)
  unint nY  : 8;    // 0..100 (8 bits)

};

struct Date2 {
  Date2(const int &pnWD, const int &pnMD,
        const int &pnM, const int &pnY):
    nWD(pnWD), nMD(pnMD), nM(pnM), nY(pnY) {}

  int nWD;
  int nMD;
  int nM;
  int nY;

};


void print(Date & d){
  std::cerr << d.nWD << " ";
  std::cerr << d.nMD << " ";
  std::cerr << d.nM << " ";
  std::cerr << d.nY << "\n";
}

void print(std::unique_ptr<Date> & d){
  std::cerr << d->nWD << " ";
  std::cerr << d->nMD << " ";
  std::cerr << d->nM << " ";
  std::cerr << d->nY << "\n";
}

int main(){

  std::vector<std::unique_ptr<Date>> a;
  std::vector<Date2> b;
  Date d(7, 0, 12, 86);
  print(d);
  std::cerr << '\n';
  a.push_back(std::make_unique<Date>(Date(7, 0, 12, 86)));
  a.push_back(std::make_unique<Date>(Date(7, 0, 12, 87)));
  b.push_back(Date2(7, 0, 12, 86));
  b.push_back(Date2(7, 0, 12, 87));
  for (auto &x: a){
    print(x);
  }
  std::cerr <<
    sizeof(std::vector<std::unique_ptr<Date>>) + (sizeof(a[0]) * a.size()) <<
    " " <<
    sizeof(std::vector<Date2>) + (sizeof(b[0]) * b.size()) <<
    '\n';

  return 0;
}
