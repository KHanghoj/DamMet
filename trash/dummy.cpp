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

void parse_vec(std::vector<std::unique_ptr<test2>> & a){
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

struct test5{
  test5(std::vector<std::shared_ptr<test2>> & _c){
    c=_c;
    std::cerr << c[0].use_count() << '\n';
  }
  std::vector<std::shared_ptr<test2>> c;
};

void parse_vec2(test4 * d){
  std::cerr << "parse_vec2" << '\n';
  for (auto &x: *(d->c)){
    std::cerr << x->a << '\n';
  }
}

void parse_vec3(test5 & d){
  std::cerr << "parse_vec3" << '\n';
  for (auto &x: d.c){
    std::cerr << x->a << '\n';
  }
}

struct deamrates{
  deamrates(std::vector<std::unique_ptr<test2>> & _d){
    d = std::move(_d);
  }
  std::vector<std::unique_ptr<test2>> d;
};

void make_a_move(std::vector<std::unique_ptr<test2>> & d){
  deamrates a(d);

  std::cerr << a.d[0]->a << '\n';
  std::cerr << "dingdong" << '\n';
  std::cerr << a.d[0]->b << '\n';
  std::cerr << d[0]->a << '\n';
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

  parse_vec(v2);
  //parse_vec(v2);

  test4 be(&v2);
  parse_vec2(&be);


  std::vector<std::shared_ptr<test2>> v3;
  v3.emplace_back(std::make_shared<test2>(1000, 10));
  v3.emplace_back(std::make_shared<test2>(10000, 100));
  std::cerr << '\n';

  std::cerr << v3[0].use_count() << '\n';

  if(1){
    test5 b2(v3);
    parse_vec3(b2);
  }
  std::cerr << v3[0].use_count() << '\n';


  make_a_move(v2);
}

  // // parsing bedfile if provided
  // if(!settings.bed_f.empty()){
  //   parse_bed_file(settings, bed_coord);
  //   settings.buffer += "\t-> Analyzing: " + std::to_string(bed_coord.size()) + " BED regions on chrom: " + settings.chrom + '\n';
  //   print_log(settings);

  //   if(bed_coord.size()==0){
  //     std::cerr << '\n' << "EXITING. NOT BED COORDS ON CHROM" << '\n';
  //     exit(EXIT_SUCCESS);
  //   }
  // }
