
// std include
#include <iostream>

// project include
#include "counter.hpp"

void usage();
int index(int argc, char* argv[]);
int count(int argc, char* argv[]);
int dump(int argc, char* argv[]);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    usage();
    
    return -1;
  }

  if (!strcmp(argv[1], "index")) {
    return index(argc, argv);
  } else if (!strcmp(argv[1], "count")) {
    return count(argc, argv);
  } else if (!strcmp(argv[1], "dump")) {
    return dump(argc, argv);
  } else {
    usage();
    return -1;
  }
  
  return -1;
  
  int k = std::stoi(argv[1]);
  int nb_threads = std::stoi(argv[5]);

  std::cerr<<"Begin of build"<<std::endl;
  set_count::Counter counter(argv[2], k, nb_threads);
  std::cerr<<"End of build"<<std::endl;

  std::cerr<<"Begin of count"<<std::endl;
  counter.count(argv[3]);
  std::cerr<<"End of count"<<std::endl;
  
  std::cerr<<"Begin of save"<<std::endl;
  counter.save(argv[4]);
  std::cerr<<"End of save"<<std::endl;
  return 0;
}

void usage() {
  std::cerr<<"USAGE:"<<std::endl;
  std::cerr<<"\tset_count [SUBCOMMAND]"<<std::endl<<std::endl;
  std::cerr<<"SUBCOMMAND:"<<std::endl;
  std::cerr<<"\tindex {kmer_size} {graph_file} {index_filename} {number_of_thread}"<<std::endl;
  std::cerr<<"\tcount {index_filename} {reads_filename} {count_filename}"<<std::endl;
  std::cerr<<"\tdump  {count_filename}"<<std::endl;
}

int index(int argc, char* argv[]) {
  if(argc != 6) {
    usage();

    return -1;
  }

  int k = std::stoi(argv[2]);
  int nb_threads = std::stoi(argv[5]);

  set_count::Counter counter(argv[3], k, nb_threads);
  
  counter.save(argv[4]);
  
  return 0;
}

int count(int argc, char* argv[]) {
  if(argc != 5) {
    usage();

    return -1;
  }

  set_count::Counter counter(argv[2]);

  counter.count(argv[3]);
  
  counter.save(argv[4]);
  
  return 0;
}

int dump(int argc, char* argv[]) {
  if(argc != 3) {
    usage();

    return -1;
  }

  set_count::Counter counter(argv[2]);

  for(auto val: *counter.count()) {
    std::cout<<set_count::kmer::kmer2seq(val.first, counter.k())<<","<<int(val.second)<<std::endl;
  }
  
  return 0;
}
