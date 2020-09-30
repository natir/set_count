#ifndef _HG_SET_COUNT_MPHF_HG_
#define _HG_SET_COUNT_MPHF_HG_

// thirdparty include
#include <BooPHF.h>
#include <kseq++/seqio.hpp>

// project inculde
#include "kmer.hpp"

namespace set_count {
  /* type definition */
  typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
  typedef boomphf::mphf<u_int64_t, hasher_t  > boophf_t;  

  /* counter object definition */
  class Mphf {
  private:
    
    std::uint8_t _k;
    kmer_t _mask;
    boophf_t* _index;
    std::vector<std::pair<std::uint64_t, std::uint8_t>> _count;

  public:

    Mphf(char* path)  {
      // open file
      std::ifstream in;
      in.open(path);

      // get kmer and init mask
      in.read(reinterpret_cast<char *>(&this->_k), sizeof(this->_k));
      this->_mask = (kmer_t(1) << (2 * this->_k)) - 1;

      // read counter value
      std::uint64_t length;
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

    Mphf(char* kmer_set, std::uint8_t k, std::uint8_t nb_threads): _k(k), _mask((kmer_t(1) << (2 * k)) - 1)  {
      // get all reference kmer
      std::vector<kmer_t> kmers;

      std::string kmer;
      std::uint64_t _count;
    
      std::ifstream in;
      in.open(kmer_set);

      while(in >> kmer >> _count) {
	kmer_t forward = kmer::seq2bit(kmer);
	kmer_t reverse = kmer::revcomp(forward, k);
      
	if (forward < reverse) {
	  kmers.push_back(forward);
	} else {
	  kmers.push_back(reverse);
	}
      }
    
      in.close();
      
      // build index
      this->_index = new boomphf::mphf<u_int64_t, hasher_t>(kmers.size(), kmers, nb_threads);

      // init counter
      this->_count = std::vector<std::pair<std::uint64_t, std::uint8_t>>(kmers.size());


      while(!kmers.empty()) {
	kmer_t kmer = kmers.back();

	u_int64_t index = this->_index->lookup(kmer);

	this->_count[index] = std::make_pair(kmer, 0);

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

      std::uint64_t length = this->_count.size();
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

#endif /* _HG_SET_COUNT_MPHF_HG_ */
