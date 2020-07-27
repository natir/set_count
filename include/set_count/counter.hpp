#ifndef _HG_SET_COUNT_COUNTER_HG_
#define _HG_SET_COUNT_COUNTER_HG_

// std include
#include <vector>
#include <fstream>
#include <algorithm>

// thirdparty include
#include <BooPHF.h>
#include <kseq++/seqio.hpp>

namespace set_count {

  /* type definition */
  typedef std::uint64_t kmer_t;
  typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
  typedef boomphf::mphf<u_int64_t, hasher_t  > boophf_t;

  /* kmer function definition */
  namespace kmer {
    kmer_t nuc2bit(char nuc) {
      return (nuc >> 1) & 0b11;
    }

    kmer_t seq2bit(std::string subseq){
      kmer_t kmer = 0;

      for(char n : subseq) {
        kmer <<= 2;
        kmer |= nuc2bit(n);
      }

      return kmer;
    }

    char bit2nuc(std::uint64_t nuc2bit) {
      switch(nuc2bit) {
      case 0: return 'A';
      case 1: return 'C';
      case 2: return 'T';
      case 3: return 'G';
      default: return 'G';
      }
    }

    std::string kmer2seq(kmer_t kmer, std::uint8_t k) {
      char buffer[31] = {0};

      for(int i = k; i > 0; i--) {
        buffer[i - 1] = bit2nuc(kmer & 0b11);

        kmer >>= 2;
      }

      std::string tmp{buffer, k};
      return tmp;
    }

    kmer_t speed_comp(kmer_t kmer) {
      return kmer ^ 0b1010101010101010101010101010101010101010101010101010101010101010;
    }

    kmer_t rev(kmer_t kmer, std::uint8_t k) {
      // Thank to needtail people ! :)
      kmer = (kmer >> 2 & 0x3333333333333333) | (kmer & 0x3333333333333333) << 2;
      kmer = (kmer >> 4 & 0x0F0F0F0F0F0F0F0F) | (kmer & 0x0F0F0F0F0F0F0F0F) << 4;
      kmer = (kmer >> 8 & 0x00FF00FF00FF00FF) | (kmer & 0x00FF00FF00FF00FF) << 8;
      kmer = (kmer >> 16 & 0x0000FFFF0000FFFF) | (kmer & 0x0000FFFF0000FFFF) << 16;
      kmer = (kmer >> 32 & 0x00000000FFFFFFFF) | (kmer & 0x00000000FFFFFFFF) << 32;
    
      return kmer >> (64 - k * 2);
    }
      
    kmer_t revcomp(kmer_t kmer, std::uint8_t k) {
      return rev(speed_comp(kmer), k);
    }
  }

  /* counter object definition */
  class Counter {
  private:

    kmer_t _mask;
    std::uint8_t _k;
    boophf_t* _index;
    std::vector<std::pair<std::uint64_t, std::uint8_t>> _count;
 
  public:

    Counter(char* path)  {
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
    
    Counter(char* kmer_set, std::uint8_t k, std::uint8_t nb_threads)  {
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

    void count(char* reads) {
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

    void save(char* path) {
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

    std::uint8_t value(std::string kmer) {
      kmer_t forward = kmer::seq2bit(kmer);
      kmer_t reverse = kmer::revcomp(forward, this->_k);

      if(forward < reverse) {
	return this->value(forward);
      } else {
	return this->value(reverse);
      }
    }
    
    std::uint8_t value(kmer_t kmer) {

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

    void inc(std::string kmer) {
      kmer_t forward = kmer::seq2bit(kmer);
      kmer_t reverse = kmer::revcomp(forward, this->_k);

      if(forward < reverse) {
	this->inc(forward);
      } else {
	this->inc(reverse);
      }
    }
    
    void inc(kmer_t kmer) {
      u_int64_t idx = this->_index->lookup(kmer);
      
      if(idx != ULLONG_MAX && this->_count[idx].first == kmer && this->_count[idx].second < UINT8_MAX) {

	this->_count[idx].second++;
      } 
    }

    std::uint8_t k()  {
      return this->_k;
    }
    
    std::vector<std::pair<std::uint64_t, std::uint8_t>>* count() {
      return &this->_count;
    }
  };
}


#endif // _HG_SET_COUNT_COUNTER_HG_