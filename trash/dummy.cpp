#include <iostream> // stdout/stdin/stderr


int main(){

  unsigned char a = 2 ;
  unsigned char b = 10 ;
  int c = 2;
  std::cerr << sizeof(a) << '\n';
  std::cerr << sizeof(c) << '\n';
  std::cerr << a + b << '\n';

  unsigned int d= 10;

  std::cerr << (a == d) << " " << (b == d) << '\n';

}
