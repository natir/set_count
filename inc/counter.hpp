#ifndef _HG_SET_COUNT_COUNTER_HG_
#define _HG_SET_COUNT_COUNTER_HG_

// std include
#include <vector>

// thirdparty include
#include <BooPHF.h>

// project include
#include "kmer.hpp"

namespace set_count {
  typedef boomphf::SingleHashFunctor<u_int64_t>  hasher_t;
  typedef boomphf::mphf<u_int64_t, hasher_t  > boophf_t;

  class Counter {
  private:

    kmer_t _mask;
    std::uint8_t _k;
    boophf_t* _index;
    std::vector<std::pair<std::uint64_t, std::uint8_t>> _count;
 
  public:

    Counter(char* path);
    Counter(char* kmer_set, std::uint8_t k, std::uint8_t nb_threads);

    ~Counter();
    
    void count(char* reads);

    void save(char* path);
    

    std::uint8_t value(std::string kmer);
    std::uint8_t value(kmer_t kmer);

    void inc(std::string kmer);
    void inc(kmer_t kmer);

    std::uint8_t k();
    std::vector<std::pair<std::uint64_t, std::uint8_t>>* count();
  };
}

#endif // _HG_SET_COUNT_COUNTER_HG_
