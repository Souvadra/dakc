#include <iostream>
#include <vector>
#include <assert.h>
#include <string>
#include <bitset>
#include <cctype>
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <map>
#include <fstream>
#include <memory>

#include <shmem.h>
#include "selector.h"

#include "common.hpp"
#include "kcounter.hpp"

#include <mpi.h>

#include <immintrin.h>

// Functions of the kmercounter class ------------------------------------------
kmer_t inline set_kmer_fast(const uint8_t *s) {
  kmer_t kmer = 0;
  for (int i = 0; i < KMERLEN; ++i) {
    kmer <<= 2;
    kmer |= static_cast<uint64_t>(s[i]);
  }
  return kmer;
}

kmer_t inline update_kmer_fast(kmer_t kmer, uint8_t s) {
  kmer <<= 2;
  kmer |= static_cast<kmer_t>(s);
  return (kmer & KMER_MASK);
}

void kmercounter::get_kmers(std::vector<kmer_t> &send_buf, 
  const uint8_t* read, int readlen, uint64_t &kmers_in_buffer) {
/*
 * Updates the kmer_send_buf with kmers extracted from the read and 
 * returns the number of kmers extracted.
 * Does NOT check for 'N' characters inside this function. Assumes 
 * reads contains only A,T,C,G characters
 */

  // define the variables 
  kmer_t curr_kmer, prv_kmer;

  // input string is smaller than a kmer 
  if (__builtin_expect(readlen < KMERLEN, 0))  return;

  // deal with the first k-mer separately (suffix only)
  curr_kmer = set_kmer_fast(read);
  send_buf[kmers_in_buffer++] = curr_kmer;
  prv_kmer = curr_kmer;

  if (__builtin_expect(readlen == KMERLEN, 0)) return;
  
  for (int i = KMERLEN; i < readlen; i++) {
    curr_kmer = update_kmer_fast(prv_kmer, read[i]);
    send_buf[kmers_in_buffer++] = curr_kmer;
    prv_kmer = curr_kmer;
  }
}

void kmercounter::read_till_buf_max(uint64_t &read_idx, 
  std::vector<kmer_t> &kmer_send_buf, bool &done_parsing, uint64_t &kmers_in_buffer) {
/*
 * Parse the read_vector and put the kmers into the kcount_buffer 
 * till either (1.) the max buffer size will exceed after adding kmers 
 * from the next read, or (2.) we exhaust the read_vector. 
 */
  char* rd = rchunk + read_idx;
  int i, left_idx;

  while (kmers_in_buffer <= (KCOUNT_BUCKET_SIZE - READLEN)) {
    // check for N characters and send the read to get_kmers function
    // then process the read and dump in the kmer_send_buf
    left_idx = 0;

    static std::vector<uint8_t> base_vec(READLEN);

    for (i = 0; i < READLEN; i++) {
      base_vec[i] = char2base(rd[i]);
      if (__builtin_expect(base_vec[i] == 0xFF, 0)) {
        get_kmers(kmer_send_buf, &base_vec[left_idx], (i - left_idx), kmers_in_buffer);
        left_idx = i + 1;
      }
    }

    if (left_idx < READLEN - 1) {
      get_kmers(kmer_send_buf, &base_vec[left_idx], (i - left_idx), kmers_in_buffer);
    }

    // move on to the next read in the (*rvec)
    read_idx += READLEN + 1; // skip current read and one \n char
    if (rchunk[read_idx] != '\0') { 
      // update the read variable
      rd = rchunk + read_idx;
    } else {
      done_parsing = true;
      break; // exits the while loop
    }
  }
}
