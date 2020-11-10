#ifndef _HG_SET_COUNT_BRISK_HG_ 
#define _HG_SET_COUNT_BRISK_HG_ 

// standard include
#include <memory>

// thirdparty include
#include <Brisk.hpp>

namespace set_count {

  class Brisk_in {
  private:
    std::uint8_t _k;
    std::uint8_t _m;
    std::unique_ptr<Brisk<std::uint8_t>> _counter;

  public:

    Brisk_in(char* path): _k(31), _m(11) {
      Parameters params(this->_k, this->_m, 4);

      this->_counter = std::make_unique<Brisk<std::uint8_t>>(params);
    }
    
    
    Brisk_in(char* kmer_set, std::uint8_t k, std::uint8_t m=11): _k(k), _m(m) {
      Parameters params(this->_k, this->_m, 4);
      this->_counter = std::make_unique<Brisk<std::uint8_t>>(params);
	 
      klibpp::KSeq record;
      klibpp::SeqStreamIn iss(kmer_set);
      std::uint8_t * data_ptr = NULL;
	
      while(iss >> record) {
	if(record.seq.length() < this->_k) {
	  continue;
	}
	
	SuperKmerEnumerator enumerator(record.seq, params.k, params.m);

	vector<kmer_full> superkmer;
	
	auto minimizer = enumerator.next(superkmer);
	while (superkmer.size() > 0) {
	  for (kmer_full & kmer : superkmer) {
	    data_ptr = this->_counter->insert(kmer);
	    *data_ptr = 0;
	  }
	  
	  // Next superkmer computation
	  superkmer.clear();
	  minimizer = enumerator.next(superkmer);
	}
      }
    }

    void count(char* reads) {
      klibpp::KSeq record;
      klibpp::SeqStreamIn iss(reads);
      uint8_t * data_ptr = NULL;
      
      while(iss >> record) {
	if(record.seq.length() < this->_k) {
	  continue;
	}

	SuperKmerEnumerator enumerator(record.seq, this->_k, this->_m);

	vector<kmer_full> superkmer;

	auto minimizer = enumerator.next(superkmer);
	while (superkmer.size() > 0) {
	  for (kmer_full & kmer : superkmer) {
	    data_ptr = this->_counter->get(kmer);
	    if(data_ptr != NULL && *data_ptr < 255) { 
	      (*data_ptr)++;
	    }
	  }
	  
	  // Next superkmer computation
	  superkmer.clear();
	  minimizer = enumerator.next(superkmer);
	}
      }
    }

    void save(char* path) {
      
    }

    std::uint8_t value(std::string kmer) {
      return this->value(str2kmer(kmer, this->_m));
    }
    
    std::uint8_t value(kmer_full kmer) {
      std::uint8_t* data_ptr = this->_counter->get(kmer);

      if(data_ptr == NULL) {
	return 0;
      }

      return *data_ptr;
    }
    
    void inc(std::string kmer) {
      this->inc(str2kmer(kmer, this->_m));
    }

    void inc(kmer_full kmer) {
      auto data_ptr = this->_counter->get(kmer);
      
      if(data_ptr != NULL && *data_ptr < 255) { 
	(*data_ptr)++;
      }
    }

    std::uint8_t k() {
      return this->_k;
    }
  };
}

#endif /* _HG_SET_COUNT_BRISK_HG_ */
