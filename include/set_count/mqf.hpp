#ifndef _HG_SET_COUNT_MQF_HG_
#define _HG_SET_COUNT_MQF_HG_

// thirdparty include
#include <kseq++/seqio.hpp>
#include <gqf.h>

// project include
#include "kmer.hpp"

namespace set_count {

  class Mqf {
  private:
    kmer_t _mask;
    std::uint8_t _k;
    QF _counter;
    std::string _path;
  
  public:
    Mqf(char* path, std::uint8_t k, double nb_uniq_kmer, char* counter_path): _k(k), _counter(), _path(counter_path) {

      /* Init MQF counter */
      std::uint64_t nslots = pow(2, ceil(log(nb_uniq_kmer * 1.05)/log(2))); 

      // nslolts is a power of two by moving bit we move the bits set at one
      // 12 is ~= to -log(10^{-5})
      qf_init(&this->_counter, nslots, 2*k, 0, 8, 0, true, counter_path, 0);
      
      /* Read uniq kmer and set in counter */
      std::string kmer;
      std::uint64_t _count;
    
      std::ifstream in;
      in.open(path);

      while(in >> kmer >> _count) {
	kmer_t forward = kmer::seq2bit(kmer);
	kmer_t reverse = kmer::revcomp(forward, k);
      
	if (forward < reverse) {
	  qf_insert(&this->_counter, forward, 1, false, false);
	} else {
	  qf_insert(&this->_counter, reverse, 1, false, false);
	}
      }
    
      in.close();
    }

    Mqf(char* path, std::uint8_t k): _k(k), _counter() {
      qf_deserialize(&this->_counter, path);
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

    void save(char* output) {
      qf_serialize(&this->_counter, output);
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
      return qf_count_key(&this->_counter, kmer);
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
      if(0 != qf_count_key(&this->_counter, kmer)) {
	qf_insert(&this->_counter, kmer, 1, false, false);
      }
    }
    
    std::uint8_t k()  {
      return this->_k;
    }
  };
}

#endif /* _HG_SET_COUNT_MQF_HG_ */
