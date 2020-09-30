#ifndef _HG_SET_COUNT_KMER_HG_
#define _HG_SET_COUNT_KMER_HG_

namespace set_count {
  typedef std::uint64_t kmer_t;

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
}

#endif /* _HG_SET_COUNT_KMER_HG_ */
