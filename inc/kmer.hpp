#ifndef _HG_SET_COUNT_KMER_HG_
#define _HG_SET_COUNT_KMER_HG_

#include <string>
#include <cstdint>

namespace set_count {

  typedef std::uint64_t kmer_t;

  namespace kmer {
    kmer_t seq2bit(std::string subseq);

    kmer_t nuc2bit(char nuc);

    std::string kmer2seq(kmer_t kmer, std::uint8_t k); 

    char bit2nuc(std::uint64_t nuc2bit);

    kmer_t canonical(kmer_t kmer, std::uint8_t k);

    bool parity_even(kmer_t kmer);
  
    kmer_t revcomp(kmer_t kmer, std::uint8_t k);

    kmer_t speed_comp(kmer_t kmer);

    kmer_t comp(kmer_t kmer, std::uint8_t k);

    kmer_t rev(kmer_t kmer, std::uint8_t k);
  }
}

#endif
