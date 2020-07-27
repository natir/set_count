// std include
#include <fstream>
#include <bitset> // debug
#include <algorithm> // debug

// thridparty
#include <kseq++/seqio.hpp>

// project include
#include "counter.hpp"

namespace set_count {

  Counter::Counter(char* path) {
    // open file
    std::ifstream in;
    in.open(path);

    // get kmer and init mask
    in.read(reinterpret_cast<char *>(&this->_k), sizeof(this->_k));
    this->_mask = (kmer_t(1) << (2 * this->_k)) - 1;

    // read counter value
    size_t length;
    in.read(reinterpret_cast<char *>(&length), sizeof(length));
    this->_count.resize(length);
    
    for(size_t i = 0; i != this->_count.size(); i++) {
      kmer_t kmer;
      std::uint8_t count;
      in.read(reinterpret_cast<char *>(&kmer), sizeof(kmer));
      in.read(reinterpret_cast<char *>(&count), sizeof(count));
      this->_count[i] = std::make_pair(kmer, count);
    }

    // read index
    this->_index = new boomphf::mphf<u_int64_t, hasher_t>();
    this->_index->load(in);

    in.close();
  }
  
  Counter::Counter(char* kmer_set, std::uint8_t k, std::uint8_t nb_threads) {
    // basic value init
    this->_k = k;

    this->_mask = (kmer_t(1) << (2 * this->_k)) - 1;

    // get all reference kmer
    std::vector<kmer_t> kmers;

    klibpp::KSeq record;
    klibpp::SeqStreamIn iss(kmer_set);
    while(iss >> record) {
      if(record.seq.length() < k) {
	continue;
      }

      kmer_t forward = kmer::seq2bit(record.seq.substr(0, k));
      kmer_t reverse = kmer::revcomp(forward, k);
      
      if(forward < reverse) {
	kmers.push_back(forward);
      } else {
	kmers.push_back(reverse);
      }

      for(char n : record.seq.substr(k)) {
	kmer_t nuc = kmer::nuc2bit(n);
	forward = ((forward << 2) & this->_mask) ^ nuc;
	reverse = (reverse >> 2) ^ ((nuc ^ 0b10)  << 2 * (k - 1));
        
	if(forward < reverse) {
	  kmers.push_back(forward);
	} else {
	  kmers.push_back(reverse);
	}
      }
    }

    // build index
    this->_index = new boomphf::mphf<u_int64_t, hasher_t>(kmers.size(), kmers, nb_threads, 1);

    // init counter
    this->_count = std::vector<std::pair<std::uint64_t, std::uint8_t>>(kmers.size());

    while(!kmers.empty()) {
      kmer_t kmer = kmers.back();

      this->_count[this->_index->lookup(kmer)] = std::make_pair(kmer, 0);

      kmers.pop_back();
    }
  }

  Counter::~Counter(){}
  
  void Counter::count(char* reads) {
    klibpp::KSeq record;
    klibpp::SeqStreamIn iss(reads);
    
    while(iss >> record) {
      if(record.seq.length() < this->_k) {
	continue;
      }
      
      kmer_t forward = kmer::seq2bit(record.seq.substr(0, this->_k));
      kmer_t reverse = kmer::revcomp(forward, this->_k);

      if(forward < reverse) {
	this->inc(forward);
      } else {
	this->inc(reverse);
      }

      for(char n : record.seq.substr(this->_k)) {
	kmer_t nuc = kmer::nuc2bit(n);
	forward = ((forward << 2) & this->_mask) ^ nuc;
	reverse = (reverse >> 2) ^ ((nuc ^ 0b10) << 2 * (this->_k - 1));
	
	if(forward < reverse) {
	  this->inc(forward);
	} else {
	  this->inc(reverse);
	}
      }
    }
  }

  
  void Counter::save(char* path) {
    std::ofstream out;
    out.open(path);

    out.write(reinterpret_cast<char const*>(&this->_k), sizeof(this->_k));

    size_t length = this->_count.size();
    out.write(reinterpret_cast<char const*>(&length), sizeof(this->_count.size()));

    for(auto val : this->_count) {
      out.write(reinterpret_cast<char const*>(&val.first), sizeof(val.first));
      out.write(reinterpret_cast<char const*>(&val.second), sizeof(val.second));
    }

    this->_index->save(out);
    
    out.close();
  }

  std::uint8_t Counter::value(std::string kmer) {
    kmer_t forward = kmer::seq2bit(kmer);
    kmer_t reverse = kmer::revcomp(forward, this->_k);

    if(forward < reverse) {
      return this->value(forward);
    } else {
      return this->value(reverse);
    }
  }
  
  std::uint8_t Counter::value(kmer_t kmer) {

    u_int64_t idx = this->_index->lookup(kmer);

    if(idx == ULLONG_MAX) {
      return 0;
    } else {
      if(this->_count[idx].first == kmer) {
	return this->_count[idx].second;
      } else {
	return 0;
      }
    }
  }


  
  void Counter::inc(std::string kmer) {
    kmer_t forward = kmer::seq2bit(kmer);
    kmer_t reverse = kmer::revcomp(forward, this->_k);

    if(forward < reverse) {
      this->inc(forward);
    } else {
      this->inc(reverse);
    }
  }
  
  void Counter::inc(kmer_t kmer) {
    u_int64_t idx = this->_index->lookup(kmer);
      
    if(idx != ULLONG_MAX && this->_count[idx].first == kmer && this->_count[idx].second < UINT8_MAX) {

      this->_count[idx].second++;
    } 
  }

  std::uint8_t Counter::k() {
    return this->_k;
  }
  
  std::vector<std::pair<std::uint64_t, std::uint8_t>>* Counter::count() {
    return &this->_count;
  }
}
