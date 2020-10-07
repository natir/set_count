#ifndef _HG_SET_COUNT_MPHF_HG_
#define _HG_SET_COUNT_MPHF_HG_

// thirdparty include
#include <BooPHF.h>
#include <kseq++/seqio.hpp>

// project inculde
#include "kmer.hpp"

#include <limits>
#include <stdexcept>
inline std::int64_t uint64_2_int64(std::uint64_t x)
{
  if (x <= std::numeric_limits<std::int64_t>::max())
    return static_cast<std::int64_t>(x);

  throw std::overflow_error("std::uint64_t value cannot be stored in a variable of type std::int64_t.");
}


namespace set_count {
  /* type definition */
  typedef boomphf::SingleHashFunctor<kmer_t>  hasher_t;
  typedef boomphf::mphf<kmer_t, hasher_t> boophf_t;  

  typedef std::uint8_t count_t;
  
  /* counter object definition */
  class Mphf {
  private:
    
    std::uint8_t _k;
    kmer_t _mask;
    boophf_t* _index;
    std::vector<kmer_t> _key;
    std::vector<count_t> _count;
    
  public:

    Mphf(char* path)  {
      // open file
      std::ifstream in;
      in.open(path);

      // get kmer and init mask
      in.read(reinterpret_cast<char *>(&this->_k), sizeof(this->_k));
      this->_mask = (kmer_t(1) << (2 * this->_k)) - 1;

      // read counter value
      std::size_t length;
      in.read(reinterpret_cast<char *>(&length), sizeof(length));

      this->_key.resize(length);
      this->_count.resize(length);

      in.read(reinterpret_cast<char*>(this->_key.data()), uint64_2_int64(length * sizeof(kmer_t)));
      in.read(reinterpret_cast<char*>(this->_count.data()), uint64_2_int64(length * sizeof(count_t)));

      // read index
      this->_index = new boomphf::mphf<u_int64_t, hasher_t>();
      this->_index->load(in);

      in.close();
    }

    Mphf(char* kmer_set, std::uint8_t k, std::uint64_t nb_uniq_kmer, std::uint8_t nb_threads): _k(k), _mask((kmer_t(1) << (2 * k)) - 1)  {
      // get all reference kmer
      std::vector<kmer_t> kmers;
      kmers.reserve(nb_uniq_kmer);
      
      std::string str_kmer;
      std::uint64_t count;
    
      std::ifstream in;
      in.open(kmer_set);

      while(in >> str_kmer >> count) {
	kmer_t forward = kmer::seq2bit(str_kmer);
	kmer_t reverse = kmer::revcomp(forward, k);
      
	if (forward < reverse) {
	  kmers.push_back(forward);
	} else {
	  kmers.push_back(reverse);
	}
      }
    
      in.close();

      // build index
      this->_index = new boomphf::mphf<u_int64_t, hasher_t>(kmers.size(), kmers, nb_threads, 1);

      // init counter
      this->_key = std::vector<kmer_t>(kmers.size());
      this->_count = std::vector<count_t>(kmers.size());
      
      while(!kmers.empty()) {
	kmer_t kmer = kmers.back();

	u_int64_t index = this->_index->lookup(kmer);

	this->_key[index] = kmer;
	this->_count[index] = count_t(0);

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
      out.write(reinterpret_cast<char const*>(&length), sizeof(size_t));
      
      out.write(reinterpret_cast<char const*>(this->_key.data()), uint64_2_int64(sizeof(kmer_t) * this->_key.size()));
      
      out.write(reinterpret_cast<char const*>(this->_count.data()), uint64_2_int64(sizeof(count_t) * this->_count.size()));
      
      this->_index->save(out);
      out.close();
    }

    count_t value(std::string kmer) {
      kmer_t forward = kmer::seq2bit(kmer);
      kmer_t reverse = kmer::revcomp(forward, this->_k);

      if(forward < reverse) {
	return this->value(forward);
      } else {
	return this->value(reverse);
      }
    }

    count_t value(kmer_t kmer) {

      u_int64_t idx = this->_index->lookup(kmer);

      if(idx == ULLONG_MAX) {
	return 0;
      } else {
	if(this->_key[idx] == kmer) {
	  return this->_count[idx];
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

      if(idx != ULLONG_MAX && this->_key[idx] == kmer && this->_count[idx] < UINT8_MAX) {
	this->_count[idx]++;
      }
    }

    std::uint8_t k()  {
      return this->_k;
    }
  };  
}

#endif /* _HG_SET_COUNT_MPHF_HG_ */
