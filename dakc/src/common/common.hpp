#ifndef __COMMON_H
#define __COMMON_H

#include <iostream>
#include <bitset>
#include <cstring>
#include <cctype>
#include <map>
#include <cinttypes>

//----------------------------
#include <shmem.h>
#include "selector.h"

#define DEBUG 0
#define PRINT 0
#define NAIVE 1

#define CURR_PE shmem_my_pe()
#define TOTAL_PE shmem_n_pes()

#define CHUNKSIZE 8 /* 512 bit registers store 8x 64 bit integers */
//-----------------------------

// Make certain variables compile time !! 
#ifndef MAX_KMER_COUNT
#define MAX_KMER_COUNT              INT64_MAX
#endif

#ifndef MIN_KMER_COUNT
#define MIN_KMER_COUNT              0
#endif 

#ifndef READLEN
#define READLEN                     150 
#endif

#ifndef KMERLEN
#define KMERLEN                     31
#endif

#ifndef HITTER
#define HITTER                      1
#endif

#ifndef BIGKSIZE
#define BIGKSIZE                    16
#endif

#ifndef KCOUNT_BUCKET_SIZE
#define KCOUNT_BUCKET_SIZE          10000
#endif

#define MINIMIZERLEN                9
#define MINCONTIGLEN                10000
// -------------------------------------

#define MAX_KMER_SIZE 64
#define KMER_T_MAX UINT64_MAX

#define KMER_MASK (((~0UL)) >> (64 - (2*KMERLEN)))
#define MINIMIZER_MASK (((~0UL)) >> (64 - (2*MINIMIZERLEN)))

typedef uint64_t kmer_t;
typedef uint64_t count_t;

typedef struct read_seq { 
    char read_data[READLEN]; 
    int read_data_size = 0;
} read_pair;

inline uint8_t char2base(const char s) {
    static const uint8_t base_map[256] = {
        /*0*/ 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
       /*16*/ 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
       /*32*/ 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
       /*48*/ 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
       /*64*/ 0xFF, 0x01, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xF0, 0xFF, 0xFF,
       /*80*/ 0xFF, 0xFF, 0xFF, 0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
       /*96*/ 0xFF, 0x01, 0xFF, 0x00, 0xFF, 0xFF, 0xFF, 0x03, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xF0, 0xFF, 0xFF,
      /*112*/ 0xFF, 0xFF, 0xFF, 0xFF, 0x02, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF,
              0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF
    };
    return base_map[static_cast<unsigned char>(s)];
}

inline char base2char(const uint8_t base_chunk) {
    static const char base_map[4] = {'C', 'A', 'T', 'G'};
    return (base_chunk < 4) ? base_map[base_chunk] : 'N';
}

inline uint8_t char2suf(const char s) {
    const char S = std::toupper(s);
    switch(S) {
        case 'A': return 0x08;
        case 'C': return 0x04;
        case 'G': return 0x02;
        case 'T': return 0x01;
    }
    // the code should never come here 
    return 0x00;
}

inline uint8_t char2pre(const char s) {
    const char S = std::toupper(s);
    switch(S) {
        case 'A': return 0x80;
        case 'C': return 0x40;
        case 'G': return 0x20;
        case 'T': return 0x10;
    }
    // the code should never come here 
    return 0x00;
}

void print_kmer(kmer_t kmer, int k); 
void print_kmer_noenter(kmer_t kmer, int k);
std::vector<bool> kmer2bitvector(kmer_t kmer, int k);

template<typename T>
void print_binary(T value) {
    const int num_bits = sizeof(T) * 8; // num of bits to store T 
    std::bitset<num_bits> binary(value);
    std::cout << binary << std::endl;
}
#endif 