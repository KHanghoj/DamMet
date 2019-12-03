#include <iostream>
#include <vector>
#include <cstdio>
#include <memory> // make_unique and uni
#include <string>
#include <cmath>  // pow
#include <iostream>
#include <fstream>
#include <sstream> //stringstream
#include <random>

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 rn_generator(rd()); //Standard mersenne_twister_engine seeded with rd()

int main(int argc, char *argv[]){

  // let me make it reproducible
  // int myseed=0;
  // rn_generator.seed(myseed);

  std::vector<std::unique_ptr<int>> my_vec;
  std::vector<std::unique_ptr<int>> my_boot;
  for(int i=0;i<10;i++)
    my_vec.emplace_back(std::make_unique<int>(i));

  // int val lrand
  std::uniform_int_distribution<int> dis(0, my_vec.size()-1);

  for(int i=0; i<10;i++){
    int idx = dis(rn_generator);
    my_boot.emplace_back(std::make_unique<int>(*my_vec[idx]));
  }

  for(const auto &val: my_boot){
    std::cerr << *val << std::endl;
  }

  return 0;
}
