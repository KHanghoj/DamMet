/* rand example: guess the number */
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <iostream> // cerr
#include <vector>

#include <chrono>
template <typename F, typename ...Args>
auto timer(F f, std::string const &label, Args && ...args) {
    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    auto holder = f(std::forward<Args>(args)...);
    auto stop = high_resolution_clock::now();
    std::cout << label << " time: " << duration_cast<nanoseconds>(stop - start).count() << "\n";
    return holder;
}
int MAX_READPOS = 20;
int STRAND=2; // 0 1
int PRIME=2;  // 0 1
int PHRED=30; // 10-40
int READPOS_MULT = STRAND * PRIME * PHRED;
int PRIME_MULT = STRAND * PHRED;
int STRAND_MULT = PHRED;

int get_idx(const int & readpos,
            const int & p,
            const int & s,
            const int & seq){
  return (readpos * READPOS_MULT +
          p * PRIME_MULT +
          s * STRAND_MULT +
          seq);
}

inline std::vector<int> idx_to_params_inline(const int & idx){
  std::div_t p = std::div(idx, READPOS_MULT);
  std::div_t pr = std::div(p.rem, PRIME_MULT);
  std::div_t strand = std::div(pr.rem, STRAND_MULT);
  return std::vector<int> {p.quot, pr.quot, strand.quot, strand.rem};
}

std::vector<int> idx_to_params(const int & idx){
  std::div_t p = std::div(idx, READPOS_MULT);
  std::div_t pr = std::div(p.rem, PRIME_MULT);
  std::div_t strand = std::div(pr.rem, STRAND_MULT);
  return std::vector<int> {p.quot, pr.quot, strand.quot, strand.rem};
}

int main ()
{

  int readpos = 4;
  int s = 1;
  int p = 1;
  int seq = 29;
  int max_size = MAX_READPOS * READPOS_MULT;
  std::vector<int> v(max_size, 0 );
  int idx = get_idx(readpos, p, s, seq);
  std::cerr << idx << " " << max_size << "\n";

  v[idx]=100;

  std::vector<int> res = idx_to_params(idx);
  for(auto &x : res){
      std::cerr << x << "\n";
  }
  std::cerr << '\n';
  std::cerr << idx << '\n';
  std::cerr << READPOS_MULT << '\n';
  std::cerr << '\n';
  int a = idx & (READPOS_MULT-1);
  std::div_t b = std::div(idx, READPOS_MULT);
  std::cerr << a << '\n';
  std::cerr << b.quot << '\n';
  std::cerr << b.rem << '\n';

  //   0 0   0  0 0 0 0 0
  // 126 64 32 16 8 4 2 1
  auto ab = timer(idx_to_params_inline,       "     inline", idx);
  auto abc = timer(idx_to_params,       "     no inline", idx);
  // https://jacksondunstan.com/articles/1946
  // Modulus
  // i % 4; // normal
  // i & 3; // bitwise [4 = 1 << 2, apply ((1 << 2) - 1), so use 3]

  return 0;
}
