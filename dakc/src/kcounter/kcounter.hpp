#ifndef __KCOUNTER_H
#define __KCOUNTER_H

#include <iostream>
#include <vector>
#include <array>
#include <unordered_map>
#include <bitset>
#include <map>

#include "common.hpp"

#define EVEN_MASK 0xAAAAAAAAAAAAAAAAULL // 101010....101010
#define ODD_MASK  0x5555555555555555ULL // 010101....010101

#define HITTERMAX 10
#define INIT_DBG_SIZE 1000000

enum MailBoxType {PUT};

enum kmer_pkt_type{NORMAL, HEAVY};

/* 
 * Since majority of the k-mers are not heavy hitters, we can send them 
 * without counting them at all. 
 * 
 * when bigk_packet.type == NORMAL, kmer_packet.count stores another k-mer 
 * when bigk_packet.type == HEAVY, kmer_packet.count stores the count of 
 * heavy hitters
 * 
 * Most of the time, the packet type will be NORMAL, so we can save a lot 
 * of space (hopefully)
 */

typedef struct kmer_packet_type {
  kmer_t kmer;
  count_t count;
} kmer_packet;

typedef struct bigk_packet_type {
  kmer_t kmers[2 * BIGKSIZE]; // second half works as 64-bit counts for heavy packets
  int size; // size is BIGKSIZE * 2 for normal, BIGKSIZE for heavy hitters
  int type;
  // uint64_t buffer; // Just to make the total packet a multiple of 64 bits !!!
} bigk_packet;

class kmer_handler: public hclib::Selector<1, bigk_packet> {
public: 
  kmer_handler(std::vector<kmer_t> *dbg, std::vector<kmer_packet> *heavydbg) 
    : dbg_(dbg), dbg_size(0), heavydbg_(heavydbg), heavydbg_size(0) {

    mb[PUT].process = [this] (bigk_packet pkt, int sender_pe) {
      this->recv_kmer(pkt, sender_pe);
    };
  }

  ~kmer_handler() {
    // Destructor code here
    dbg_->resize(dbg_size);
    #if HITTER
    heavydbg_->resize(heavydbg_size);
    #endif
  }

private: 
  std::vector<kmer_t> *dbg_;
  std::vector<kmer_packet> *heavydbg_;
  uint32_t dbg_size, heavydbg_size;
  void recv_kmer(bigk_packet pkt, int sender_pe);
};

// kmer counting class
class kmercounter {
private:
public:
  std::vector<kmer_t> *vectordbg;
  std::vector<kmer_packet> *heavydbg;
  std::vector<kmer_packet> *lightdbg;
  char* rchunk;

  const uint8_t pre_delete_mask[4] = {0x7F, 0xBF, 0xDF, 0xEF};
  const uint8_t suf_delete_mask[4] = {0xF7, 0xFB, 0xFD, 0xFE};

  kmercounter(char* read_chunk, std::vector<kmer_t> &vectordbg) {
    
    this->rchunk = read_chunk;
    this->vectordbg = &vectordbg;
    this->vectordbg->resize(INIT_DBG_SIZE);

    #if HITTER
    this->heavydbg = new std::vector<kmer_packet>();
    this->heavydbg->resize(INIT_DBG_SIZE);
    #endif

    this->lightdbg = new std::vector<kmer_packet>();
    this->lightdbg->resize(INIT_DBG_SIZE);

    perform_kcount();
  }

  ~kmercounter() {
    #if HITTER
    delete heavydbg;
    #endif
  }

  void get_kmers(std::vector<kmer_t> &sendbuf, const uint8_t* read, int readlen, uint64_t &kmers_in_buffer);
  void read_till_buf_max(uint64_t &read_idx, std::vector<kmer_t> &sendbuf, bool &done_parsing, uint64_t &kmers_in_buffer);
  void flush_buffer(std::vector<kmer_t> &kcount_buffer, uint64_t &kmers_in_buffer, kmer_handler* kmer_selector, 
    std::vector<bigk_packet> &heavy_send_pkt_vec, std::vector<bigk_packet> &big_send_pkt_vec);

  void perform_kcount();
};

#endif 
