#include <string>
#include <vector>
#include <iostream>


int main(int argc, char *argv[]){
  std::vector<std::string> args;
  for(int i=1;i<argc;i++){
    args.push_back(std::string(argv[i]));
  }

  for(size_t i=0;i<args.size();i++){
    if(args[i] == "-b" || args[i].compare("-bam")==0){
      i++;
      std::cerr << args[i] << std::endl;
    }
    if(args[i] == "-r" || args[i].compare("-ref")==0){
      i++;
      std::cerr << args[i] << std::endl;
    }
  }

  std::cerr << "break\n";

  for(auto i=args.begin();i!=args.end();i++){
    if(*i == "-bam" ){
      i++;
      std::cerr << *i << std::endl;
    }
    if((*i) == "-ref"){
      i++;
      std::cerr << (*i) << std::endl;
    }
  }

  return 0;
}
