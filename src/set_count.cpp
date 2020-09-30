
// std include
#include <iostream>

// project include
#include "set_count.hpp"

void usage();
int mphf_index(int argc, char* argv[]);
int mphf_count(int argc, char* argv[]);
int mphf_dump(int argc, char* argv[]);
int mqf_index(int argc, char* argv[]);
int mqf_count(int argc, char* argv[]);
int mqf_dump(int argc, char* argv[]);

int main(int argc, char* argv[]) {
  if (argc < 2) {
    usage();

    return -1;
  }

  if (!strcmp(argv[1], "mphf_index")) {
    return mphf_index(argc, argv);
  } else if (!strcmp(argv[1], "mphf_count")) {
    return mphf_count(argc, argv);
  } else if (!strcmp(argv[1], "mphf_dump")) {
    return mphf_dump(argc, argv);
  } else if (!strcmp(argv[1], "mqf_index")) {
    return mqf_index(argc, argv);
  } else if (!strcmp(argv[1], "mqf_count")) {
    return mqf_count(argc, argv);
  } else if (!strcmp(argv[1], "mqf_dump")) {
    return mqf_dump(argc, argv);
  } else {
    usage();
    return -1;
  }

  usage();
  return -1;
}

void usage() {
  std::cerr<<"USAGE:"<<std::endl;
  std::cerr<<"\tset_count [SUBCOMMAND]"<<std::endl<<std::endl;
  std::cerr<<"SUBCOMMAND:"<<std::endl;
  std::cerr<<"\tmphf_index {kmer_size} {uniq_kmer_file} {index_filename} {number_of_thread}"<<std::endl;
  std::cerr<<"\tmphf_count {index_filename} {reads_filename} {count_filename}"<<std::endl;
  std::cerr<<"\tmphf_dump  {count_filename} {reference_filename}"<<std::endl;
  std::cerr<<"\tmqf_index {kmer_size} {number_of_uniq_kmer} {uniq_kmer_file} {mqf_save}"<<std::endl;
  std::cerr<<"\tmqf_count {kmer_size} {reads_filename} {count_filename}"<<std::endl;
  std::cerr<<"\tmqf_dump {kmer_size} {reads_filename} {reference_filename}"<<std::endl;
}

int mphf_index(int argc, char* argv[]) {
  if(argc != 6) {
    usage();

    return -1;
  }

  std::uint8_t k = std::uint8_t(std::stoi(argv[2]));

  std::uint8_t nb_threads = std::uint8_t(std::stoi(argv[5]));

  set_count::Mphf counter(argv[3], k, nb_threads);

  counter.save(argv[4]);

  return 0;
}

int mphf_count(int argc, char* argv[]) {
  if(argc != 5) {
    usage();

    return -1;
  }

  set_count::Mphf counter(argv[2]);

  counter.count(argv[3]);

  counter.save(argv[4]);

  return 0;
}

int mphf_dump(int argc, char* argv[]) {
  if(argc != 4) {
    usage();

    return -1;
  }

  set_count::Mphf counter(argv[2]);
  set_count::kmer_t mask = (set_count::kmer_t(1) << (2 * counter.k())) - 1;
  
  klibpp::KSeq record;
  klibpp::SeqStreamIn iss(argv[3]);

  while(iss >> record) {
    if(record.seq.length() < counter.k()) {
      continue;
    }

    set_count::kmer_t forward = set_count::kmer::seq2bit(record.seq.substr(0, counter.k()));
    set_count::kmer_t reverse = set_count::kmer::revcomp(forward, counter.k());

    if(forward < reverse) {
	counter.value(forward);
      } else {
	counter.value(reverse);
    }

    for(char n : record.seq.substr(counter.k())) {
      set_count::kmer_t nuc = set_count::kmer::nuc2bit(n);
      forward = ((forward << 2) & mask) ^ nuc;
      reverse = (reverse >> 2) ^ ((nuc ^ 0b10)  << 2 * (counter.k() - 1));

      if(forward < reverse) {
	std::cout<<counter.value(forward)<<std::endl;
      } else {
	std::cout<<counter.value(reverse)<<std::endl;
      }
    }
  }
	
  return 0;
}

int mqf_index(int argc, char* argv[]) {
  if(argc != 6) {
    usage();

    return -1;
  }
  
  std::uint8_t k = std::uint8_t(std::stoi(argv[2]));
  double nb_uniq_kmer = std::stod(argv[3]);

  set_count::Mqf counter(argv[4], k, nb_uniq_kmer, argv[5]);

  counter.save(argv[5]);

  return 0;
}

int mqf_count(int argc, char* argv[]) {
  if(argc != 5) {
    usage();

    return -1;
  }
  
  std::uint8_t k = std::uint8_t(std::stoi(argv[2]));
  
  set_count::Mqf counter(argv[4], k);

  counter.count(argv[3]);

  counter.save(argv[4]);

  return 0;
}

int mqf_dump(int argc, char* argv[]) {
  if(argc != 5) {
    usage();

    return -1;
  }

  std::uint8_t k = std::uint8_t(std::stoi(argv[2]));
  set_count::Mqf counter(argv[3], k);
  set_count::kmer_t mask = (set_count::kmer_t(1) << (2 * counter.k())) - 1;
  
  klibpp::KSeq record;
  klibpp::SeqStreamIn iss(argv[4]);

  while(iss >> record) {
    if(record.seq.length() < counter.k()) {
      continue;
    }

    set_count::kmer_t forward = set_count::kmer::seq2bit(record.seq.substr(0, counter.k()));
    set_count::kmer_t reverse = set_count::kmer::revcomp(forward, counter.k());

    if(forward < reverse) {
	counter.value(forward);
      } else {
	counter.value(reverse);
    }

    for(char n : record.seq.substr(counter.k())) {
      set_count::kmer_t nuc = set_count::kmer::nuc2bit(n);
      forward = ((forward << 2) & mask) ^ nuc;
      reverse = (reverse >> 2) ^ ((nuc ^ 0b10)  << 2 * (counter.k() - 1));

      if(forward < reverse) {
	std::cout<<counter.value(forward)<<std::endl;
      } else {
	std::cout<<counter.value(reverse)<<std::endl;
      }
    }
  }
	
  return 0;
}
