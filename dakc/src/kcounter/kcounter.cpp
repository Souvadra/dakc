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

#ifdef PROFILE
#define ENABLE_TRACE
#endif

#include "selector.h"

#include "common.hpp"
#include "kcounter.hpp"
#include "ska_sort.hpp"

#include <mpi.h>
#include <immintrin.h>

uint64_t MurmurHash64A (uint64_t key, uint64_t seed) {
  const uint64_t m = 0xc6a4a7935bd1e995;
  const int r = 47;

  uint64_t h = seed ^ (8 * m);

  uint64_t k = key;
  k *= m; 
  k ^= k >> r; 
  k *= m; 
  
  h ^= k;
  h *= m; 

  h ^= h >> r; 
  h *= m; 
  h ^= h >> r; 

  return h;
}

inline int owner_pe(kmer_t kmer) {
  /* example of a randomly chosen 64-bit seed */
  const uint64_t seed = 0x9E3779B97F4A7C15;
  return MurmurHash64A(kmer, seed) % TOTAL_PE;
}

int binary_search(const std::vector<kmer_packet> &arr, kmer_t kmer, int &left_, int right) {

  auto it = std::lower_bound(arr.begin() + left_, arr.begin() + right + 1, kmer, 
    [](const kmer_packet &pkt, kmer_t kmer) {
      return pkt.kmer < kmer;
  });

  if (it != arr.begin() + right + 1 && it->kmer == kmer) {
    left_ = std::distance(arr.begin(), it);
    return left_;
  }

  return -1;
}

void sort_and_merge_duplicate_kmer_packets(std::vector<kmer_packet> &vec, uint32_t &size) {

  if (__builtin_expect(size == 0, 0)) return;

  /* First sort the vector */
  ska_sort(vec.begin(), vec.begin() + size, [](const kmer_packet &a) {return a.kmer;});

  /* using iterators to make the code more efficient */
  auto it = vec.begin();
  auto end = vec.begin() + size;

  /* output iterator, location at which we'll write */
  auto out_it = vec.begin();

  kmer_packet curr_pkt = std::move(*it);
  it++;

  while (it != end) {
    if (it->kmer == curr_pkt.kmer) {
      curr_pkt.count += it->count;
    } else {
      /* write it to the output vector (same as input vector) */
      *out_it = std::move(curr_pkt);
      out_it++;

      /* update the curr_pkt variable */
      curr_pkt = std::move(*it);
    }
    it++;
  }

  *out_it = std::move(curr_pkt);
  size = std::distance(vec.begin(), out_it) + 1;
}

#if 0
void simd_transfer(const kmer_t* src, kmer_t* dest, 
  uint32_t dest_start, uint32_t src_size) {
  
  uint32_t i = 0;
  for (; i < src_size; i += CHUNKSIZE) {
    __m512i src_chunk = _mm512_loadu_si512(src + i);
    _mm512_storeu_si512(dest + dest_start + i, src_chunk);
  }

  for (; i < src_size; i++) {
    dest[dest_start + i] = src[i];
  }
}
#endif

// Message Handler -------------------------------------------------------------
void kmer_handler::recv_kmer(bigk_packet pkt, int sender_pe) {
  if (__builtin_expect(pkt.type == NORMAL, 1)) {
    if (__builtin_expect(dbg_size + pkt.size > dbg_->size(), 0)) {
      dbg_->resize(2 * dbg_size);
    }
    // simd_transfer(&pkt.kmers[0], dbg_->data(), dbg_size, pkt.size);
    for (int i = 0; i < pkt.size; i++) {
      (*dbg_)[dbg_size + i] = pkt.kmers[i];
    }
    dbg_size += pkt.size;
  } else { // HEAVY HITTER TYPE PACKET
    if (__builtin_expect(heavydbg_size + pkt.size > heavydbg_->size(), 0)) {
      heavydbg_->resize(2 * heavydbg_size);
    }

    for (int i = 0; i < pkt.size; i++) {
      (*heavydbg_)[heavydbg_size + i] = {pkt.kmers[i], pkt.kmers[BIGKSIZE + i]};
    }
    
    heavydbg_size += pkt.size;
  }
}

void init_packets(std::vector<bigk_packet> &pkt_vec, int type) {
  for (int i = 0; i < TOTAL_PE; i++) {
    pkt_vec[i].size = 0;
    pkt_vec[i].type = type;
  }
}

void empty_packets(std::vector<bigk_packet> &pkt_vec, kmer_handler* kmer_selector) {
  for (int i = 0; i < TOTAL_PE; i++) {
    if (pkt_vec[i].size > 0) {
      kmer_selector->send(PUT, pkt_vec[i], i);
    }
    pkt_vec[i].size = 0;
  }
}

void inline add_in_normal_packet(std::vector<bigk_packet_type> &normal_vec, kmer_t kmer, 
    kmer_handler* kmer_selector) {
  int owner = owner_pe(kmer);
  bigk_packet &bigpkt = normal_vec[owner];

  bigpkt.kmers[bigpkt.size] = kmer;
  bigpkt.size++;

  if (bigpkt.size == BIGKSIZE * 2) {
    kmer_selector->send(PUT, bigpkt, owner);
    bigpkt.size = 0;
  }
}

void inline add_in_heavy_packet(std::vector<bigk_packet_type> &heavy_vec, kmer_t kmer, 
    count_t count, kmer_handler* kmer_selector) {
  int owner = owner_pe(kmer);
  bigk_packet &bigpkt = heavy_vec[owner];

  bigpkt.kmers[bigpkt.size] = kmer;
  bigpkt.kmers[BIGKSIZE + bigpkt.size] = count;
  bigpkt.size++;

  if (bigpkt.size == BIGKSIZE) {
    kmer_selector->send(PUT, bigpkt, owner);
    bigpkt.size = 0;
  }
}

void send2sendbuf(kmer_t curr_kmer, count_t curr_count, kmer_handler* kmer_selector, 
  std::vector<bigk_packet_type> &hitter_vec, std::vector<bigk_packet_type> &normal_vec) {

  #if DEBUG 
  assert(curr_count > 0);
  #endif
  
  switch (curr_count) {
    case 1:
      add_in_normal_packet(normal_vec, curr_kmer, kmer_selector);
      break;
    case 2:
      add_in_normal_packet(normal_vec, curr_kmer, kmer_selector);
      add_in_normal_packet(normal_vec, curr_kmer, kmer_selector);
      break;
    default:
      add_in_heavy_packet(hitter_vec, curr_kmer, curr_count, kmer_selector);
      break;
  }
}

void kmercounter::flush_buffer(std::vector<kmer_t> &kcount_buffer, uint64_t &kmers_in_buffer,
    kmer_handler* kmer_selector,std::vector<bigk_packet_type> &hitter_vec, 
    std::vector<bigk_packet_type> &normal_vec) {
  
  int i, owner; 
  kmer_t kmer;

  #if !HITTER
  for (i = 0; i < kmers_in_buffer; i++) {
    kmer = kcount_buffer[i];
    add_in_normal_packet(normal_vec, kmer, kmer_selector);
  }
  #endif

  #if HITTER
  ska_sort(kcount_buffer.begin(), kcount_buffer.begin() + kmers_in_buffer);
  kmer_t curr_kmer = kcount_buffer[0];
  count_t curr_count = 1;

  /* For every packet I send to the heavy buffer, I send a single packet to the normal buffer */
  /* Hence, all the kmers in heavy dbg should already be there in the normal dbg*/
  for (i = 1; i < kmers_in_buffer; i++) {
    if (kcount_buffer[i] == curr_kmer) {
      curr_count++;
    } else {
      send2sendbuf(curr_kmer, curr_count, kmer_selector, hitter_vec, normal_vec);
      curr_kmer = kcount_buffer[i];
      curr_count = 1;
    }
  }
  send2sendbuf(curr_kmer, curr_count, kmer_selector, hitter_vec, normal_vec);
  
  #endif
}

void kmercounter::perform_kcount() {
/*
 * the main function of kmercounter class that takes the input vector 
 * and build the de bruijn graph in terms of a lookup table (implicitly)
 */ 
  double starttime, endtime, localtime, globaltime;

  if (CURR_PE == 0) {
    #if HITTER
    std::cout << "HITTER flag in ON " << std::endl;
    #else 
    std::cout << "HITTER flag in OFF " << std::endl;
    #endif
  }

  starttime = MPI_Wtime();
  kmer_handler* kmer_selector = new kmer_handler(vectordbg, heavydbg);

  hclib::finish([=]() {
    bool done_parsing = false;
    // initialize the variables
    uint64_t kmers_in_buffer = 0;
    std::vector<kmer_t> kcount_buffer(KCOUNT_BUCKET_SIZE + (2 * READLEN));
    std::vector<bigk_packet> big_send_pkt_vec(TOTAL_PE);
    
    std::vector<bigk_packet> heavy_send_pkt_vec;
    #if HITTER 
    heavy_send_pkt_vec.resize(TOTAL_PE);
    #endif
    
    uint64_t read_idx = 0;
    uint64_t iteration = 0;

    init_packets(big_send_pkt_vec, NORMAL);

    #if HITTER
    init_packets(heavy_send_pkt_vec, HEAVY);
    #endif

    // start the kmer parsing and sending to its owner process
    kmer_selector->start();
    while (!done_parsing) {
      read_till_buf_max(read_idx, kcount_buffer, done_parsing, kmers_in_buffer);
      flush_buffer(kcount_buffer, kmers_in_buffer, kmer_selector, heavy_send_pkt_vec, big_send_pkt_vec);
      kmers_in_buffer = 0;
    }
    empty_packets(big_send_pkt_vec, kmer_selector);

    #if HITTER 
    empty_packets(heavy_send_pkt_vec, kmer_selector);
    #endif 

    kmer_selector->done(PUT);
  });

  #ifdef ENABLE_TRACE
  // std::cout << "PE: " << hclib::TOTAL_L3_MISSES_SOUVI << " | L3 misses" << std::endl;

  long long localL3Misses = hclib::TOTAL_L3_MISSES_SOUVI;
  long long globalL3Misses = 0;

  MPI_Reduce(&localL3Misses, &globalL3Misses, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  globalL3Misses = globalL3Misses * 24 / TOTAL_PE;

  if (CURR_PE == 0) {
    std::cout << "Phase 1 L3 misses: " << globalL3Misses << std::endl;
  }
  #endif

  #ifdef ENABLE_TCOMM_PROFILING
  char profile_name[] = "kmer_counting";
  kmer_selector->print_profiling(profile_name);
  #endif
  delete kmer_selector;

  uint32_t low_freq_size = 0;
  uint32_t vectordbg_size = vectordbg->size();

  #if HITTER
  uint32_t high_freq_size = heavydbg->size();

  /* First, sort and merge the duplicates in the high frequency arrays */
  sort_and_merge_duplicate_kmer_packets(*heavydbg, high_freq_size);

  /* Now, deal with the low frequency kmer array */
  ska_sort(vectordbg->begin(), vectordbg->begin() + vectordbg_size);

  kmer_t curr_kmer = (*vectordbg)[0];
  count_t curr_count = 1;
  int idx;
  int left = 0, right = high_freq_size - 1;

  #ifdef BENCHMARK
  uint32_t binary_search_hit = 0;
  #endif

  for (uint32_t i = 1; i < vectordbg_size; i++) {
    kmer_t kmer = (*vectordbg)[i];
    if (kmer == curr_kmer) {
      curr_count++;
    } else {
      idx = binary_search(*heavydbg, kmer, left, right);
      // int idx2 = binary_search_vanilla(*heavydbg, kmer, high_freq_size);

      // if (idx != idx2) {
      //   std::cout << "left: " << left << " right: " << right  << " idx: " << idx << " idx2: " << idx2 << " kmer: " << kmer << std::endl;
      // }
      
      if (__builtin_expect(idx != -1, 0)) {
        (*heavydbg)[idx].count += curr_count;
        #ifdef BENCHMARK_
        binary_search_hit++;
        #endif
      } else {
        (*lightdbg)[low_freq_size] = {curr_kmer, curr_count};
        low_freq_size++;

        if (__builtin_expect(low_freq_size == (*lightdbg).size(), 0)) {
          (*lightdbg).resize(2 * low_freq_size);
        }

      }
      curr_kmer = (*vectordbg)[i];
      curr_count = 1;
    }
  }

  /* Now, (*lightdbg) and heavydbg are two sorted arrays that contain all the k-mers 
    and their counts in the input dataset. For querying, give priority to the heavy- 
    dbg first, and if a k-mer is not present there, then look for that in (*lightdbg) */

  #else // HITTER == 0
  
  /* Sort and merge the duplicates in the normal kmer array */
  ska_sort(vectordbg->begin(), vectordbg->begin() + vectordbg_size);
  kmer_t curr_kmer = (*vectordbg)[0];
  count_t curr_count = 1;

  for (uint32_t i = 1; i < vectordbg_size; i++) {
    kmer_t kmer = (*vectordbg)[i];
    if (kmer == curr_kmer) {
      curr_count++;
    } else {
      (*lightdbg)[low_freq_size] = {curr_kmer, curr_count};
      low_freq_size++;

      if (__builtin_expect(low_freq_size == (*lightdbg).size(), 0)) {
        (*lightdbg).resize(2 * low_freq_size);
      }

      curr_kmer = (*vectordbg)[i];
      curr_count = 1;
    }
  }

  /* Now, just query sorted (*lightdbg) array to get all the k-mers and their counts */
  #endif

  endtime = MPI_Wtime();

  vectordbg->clear(); // free the memory

  localtime = endtime - starttime; 
  MPI_Reduce(&localtime, &globaltime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  #if DEBUG 
  std::cout << "PE: " << CURR_PE << " | final_kmer_list_size: " << low_freq_size << std::endl;
  #endif

  if (CURR_PE == 0){
    std::cout << "kmer counting time: " << globaltime << " seconds"; 
    std::cout << std::endl;
  }

  #ifdef BENCHMARK
  // #if 1
    uint64_t local_distinct_kmers = low_freq_size;
    uint64_t local_kmers = 0;

    uint64_t global_kmers, global_distinct_kmers;

    for (int i = 0; i < low_freq_size; i++) {
      local_kmers += (*lightdbg)[i].count;
    }

    #if HITTER 
    local_distinct_kmers += high_freq_size;
    for (int i = 0; i < high_freq_size; i++) {
      local_kmers += (*heavydbg)[i].count;
    }

    uint64_t lnormal_size, lheavy_size, gnormal_size, gheavy_size;
    lnormal_size = low_freq_size;
    lheavy_size = high_freq_size;

    MPI_Reduce(&lnormal_size, &gnormal_size, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&lheavy_size, &gheavy_size, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);

    uint64_t total_heavy_size = 0;
    MPI_Reduce(&lheavy_size, &total_heavy_size, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    uint32_t global_binary_search_hits;
    MPI_Reduce(&binary_search_hit, &global_binary_search_hits, 1, MPI_UINT32_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (CURR_PE == 0) {
      std::cout << "max normal size: " << gnormal_size << std::endl;
      std::cout << "max heavy size: " << gheavy_size << std::endl;
      std::cout << "total heavy size: " << total_heavy_size << std::endl;
      std::cout << "global binary search hits: " << global_binary_search_hits << std::endl;
    }
    #endif

    // std::cout << "PE: " << CURR_PE << " Local kmers: " << local_kmers << std::endl;
    // std::cout << "PE: " << CURR_PE << " Local distinct kmers: " << local_distinct_kmers << std::endl;

    MPI_Allreduce(&local_distinct_kmers, &global_distinct_kmers, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_kmers, &global_kmers, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

    uint64_t global_max_kmers, global_max_distinct_kmers;
    MPI_Allreduce(&local_distinct_kmers, &global_max_distinct_kmers, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&local_kmers, &global_max_kmers, 1, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);

    double mean_dist = (double)global_distinct_kmers / (double)TOTAL_PE;
    double mean  = (double)global_kmers / (double)TOTAL_PE;

    double local_variable = (local_kmers - mean) * (local_kmers - mean);
    double local_dist_variable = (local_distinct_kmers - mean_dist) * (local_distinct_kmers - mean_dist);

    double variance, dist_variance;

    MPI_Reduce(&local_variable, &variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_dist_variable, &dist_variance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    variance /= (double)TOTAL_PE;
    dist_variance /= (double)TOTAL_PE;

    variance = sqrt(variance);
    dist_variance = sqrt(dist_variance);

    if (CURR_PE == 0) {
      std::cout << "global max kmers: " << global_max_kmers << std::endl;
      std::cout << "global max distinct kmers: " << global_max_distinct_kmers << std::endl;
      std::cout << "global_kmers: " << global_kmers << std::endl;
      std::cout << "global_distinct_kmers: " << global_distinct_kmers << std::endl;
      std::cout << "std_dev: " << variance << std::endl;
      std::cout << "dist_std_dev: " << dist_variance << std::endl;
    }

  #endif
}
