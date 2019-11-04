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
template <typename F, typename ...Args>
auto timer2(F f, std::string const &label, Args && ...args) {
    using namespace std::chrono;

    auto start = high_resolution_clock::now();
    f(std::forward<Args>(args)...);
    auto stop = high_resolution_clock::now();
    std::cout << label << " time: " << duration_cast<nanoseconds>(stop - start).count() << "\n";
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

void my_divmod(const int & a, const int &b, int & res, int &rem){
  res = a/b;
  rem = a%b;
}

int my_divmod(const int & a, const int &b, int & res){
  res = a/b;
  return(a%b);
}


void idx_to_params_kh2(const int & idx, std::vector<int> & res){
  int rem, rem1;
  my_divmod(idx, READPOS_MULT, res[0], rem);
  my_divmod(rem, PRIME_MULT, res[1], rem1);
  my_divmod(rem1, STRAND_MULT, res[2], res[3]);
}
void idx_to_params_kh3(const int & idx, std::vector<int> & res){
  int rem, rem1;
  rem = my_divmod(idx, READPOS_MULT, res[0]);
  rem = my_divmod(rem, PRIME_MULT, res[1]);
  res[3] = my_divmod(rem, STRAND_MULT, res[2]);
}

std::vector<int> idx_to_params_kh(const int & idx){
  int rem, rem1;
  std::vector<int> res(4,0);
  my_divmod(idx, READPOS_MULT, res[0], rem);
  my_divmod(rem, PRIME_MULT, res[1], rem1);
  my_divmod(rem1, STRAND_MULT, res[2], res[3]);
  return(res);
}

int main ()
{

  int readpos = 4;
  int p = 1;
  int s = 0;
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
  int a = idx & (READPOS_MULT);
  std::div_t b = std::div(idx, READPOS_MULT);
  std::cerr << a << '\n';
  std::cerr << b.quot << '\n';
  std::cerr << b.rem << '\n';

  // //   0 0   0  0 0 0 0 0
  // // 126 64 32 16 8 4 2 1
  std::vector<int> res2(4, 0);
  auto ab = timer(idx_to_params_inline,       "     inline", idx);
  auto abc = timer(idx_to_params,       "     no inline", idx);
  auto abcd = timer(idx_to_params_kh,       "     kh", idx);
  timer2(idx_to_params_kh2,       "     kh2", idx, res2);
  timer2(idx_to_params_kh3,       "     kh3", idx, res2);
  std::vector<int> res3 = idx_to_params_kh(idx);
  for (int & x1 : res3){
    std::cerr << x1 << "\n";
  }


  idx_to_params_kh2(idx, res2);
  for (int & x1 : res2){
    std::cerr << x1 << "\n";
  }

  idx_to_params_kh3(idx, res2);
  for (int & x1 : res2){
    std::cerr << x1 << "\n";
  }

  // https://jacksondunstan.com/articles/1946
  // Modulus
  // i % 4; // normal
  // i & 3; // bitwise [4 = 1 << 2, apply ((1 << 2) - 1), so use 3]

  return 0;
}
