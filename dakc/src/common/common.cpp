#include <iostream>
#include <limits>
#include <cstdint>

#include <shmem.h>
#include "selector.h"

#include "common.hpp"

#if DEBUG
void print_kmer(kmer_t kmer, int k) { 
  int left_zeros = MAX_KMER_SIZE - (2 * k);
  kmer <<= left_zeros;
  for (int i = 0; i < k; i++) { 
      uint8_t base = static_cast<uint8_t>((kmer >> (MAX_KMER_SIZE - 2)));
      std::cout << base2char(base);
      kmer <<= 2;
  }
  std::cout << std::endl;
}

void print_kmer_noenter(kmer_t kmer, int k) {
  int left_zeros = MAX_KMER_SIZE - (2 * k);
  kmer <<= left_zeros;
  for (int i = 0; i < k; i++) { 
    uint8_t base = static_cast<uint8_t>((kmer >> (MAX_KMER_SIZE - 2)));
    std::cout << base2char(base);
    kmer <<= 2;
  }
}

std::vector<bool> kmer2bitvector(kmer_t kmer, int k) {
  int left_zeros = MAX_KMER_SIZE - (2 * k);
  kmer <<= left_zeros;

  std::vector<bool> output(2 * k);
  for (int i = 0; i < 2*k; i++) {
    output[i] = static_cast<bool>((kmer >> 63) & 1);
    kmer <<= 1;
  }

  assert(kmer == 0);

  return output;
}
#endif
