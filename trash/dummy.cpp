#include <iostream> // stdout/stdin/stderr
#include <memory>
#include <vector>
struct test{
  test(const int & _a,
       std::unique_ptr<int> & _b,
       std::vector<std::unique_ptr<int>> &_c):
    a(_a), b(_b) {
    c=std::move(_c);
  }
  int a;
  std::unique_ptr<int> &b;
  std::vector<std::unique_ptr<int>> c;
};


struct test2{
  test2(const int&_a,
        const int&_b):
    a(_a), b(_b) {}
  int a,b;
};

void test3(std::vector<std::unique_ptr<test2>> & a){
  std::cerr << "test" << '\n';
  for (auto&x: a){
    std::cerr << x->a << '\n';
  }
}

struct test4{
  test4(std::vector<std::unique_ptr<test2>> * _c){
    c=_c;
  }
  std::vector<std::unique_ptr<test2>> * c;
};

void test5(test4 * d){
  std::cerr << "test5" << '\n';
  for (auto &x: *(d->c)){
    std::cerr << x->a << '\n';
  }

}

int main(){

  unsigned char a = 2 ;
  unsigned char b = 10 ;
  int c = 2;
  std::cerr << sizeof(a) << '\n';
  std::cerr << sizeof(c) << '\n';
  std::cerr << a + b << '\n';

  unsigned int d= 10;

  std::cerr << (a == d) << " " << (b == d) << '\n';

  std::unique_ptr<int>  p1 = std::make_unique<int>(100);
  std::vector<std::unique_ptr<int>> v1;
  v1.emplace_back(std::make_unique<int>(1000));

  std::cerr << *v1[0] << '\n';
  test s(1, p1, v1); // , v1);
  for (std::unique_ptr<int> & x : s.c){
    std::cerr << s.a + *s.b + *x << '\n';
  }

  std::cerr << s.a + *s.b << '\n';

  std::vector<std::unique_ptr<test2>> v2;
  v2.emplace_back(std::make_unique<test2>(1000, 10));
  v2.emplace_back(std::make_unique<test2>(10000, 100));

  for (const std::unique_ptr<test2> &x: v2){
    std::cerr << (*x).a << '\n';
    std::cerr << x->b << '\n';
  }

  test3(v2);
  test3(v2);

  test4 be(&v2);
  test5(&be);

}
