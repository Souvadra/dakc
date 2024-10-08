// ***********************************************************************
//
// PaKman: Algorithm for generating genomic contigs on distributed-memory machines
// 
// Priyanka Ghosh (Pacific Northwest National Laboratory)
// Sriram Krishnamoorthy (Pacific Northwest National Laboratory)
// Ananth Kalyanaraman (Washington State University)
//               
//
// ***********************************************************************
//
//       Copyright (2020) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************


/* 
 * Future Editor: Souvadra Hati  
 * Comments: Removed the unnecessary assert statements and included ska_sort 
 *           as a replacement for std::sort. Improved performance by 2x.
 */

#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
#include "timers.h"
#include "ska_sort.hpp"

extern long int MAX_KMER_COUNT;
extern int rank, size;
extern int coverage;
extern int num_threads;
extern int num_batch_transfers;
uint64_t num_recalculate_lmer, global_num_recalculate_lmer;

extern std::vector<lmer_t> lmer_frequency;
extern std::vector<lmer_t> global_lmer_frequency;
extern std::vector<KmerPairs> kmer_proc_buf;

void sort_recv_buffer(
    std::vector<KmerPairs>& kmer_recv_buf, 
    std::vector<int>& rcounts_kmer, 
    std::vector<int>& rdisp_kmer
) {
    size_t rsize=0;
    std::vector<int> offsets(rdisp_kmer);
    offsets.push_back(rcounts_kmer[size-1]+rdisp_kmer[size-1]);
    for (int t=0; t<rcounts_kmer.size(); t++) rsize += rcounts_kmer[t];

    while(offsets.size()>2) {
        std::vector<int> new_offsets;
        int x = 0;
        while(x+2 < offsets.size()) {
            std::inplace_merge(
                kmer_recv_buf.begin()+offsets[x], kmer_recv_buf.begin()+offsets[x+1], 
                kmer_recv_buf.begin()+offsets[x+2],  // this *might* be at the end
                [](const auto& i, const auto& j) {return i.seq < j.seq;} 
            );

            /* now they are sorted, we just put offsets[x] and offsets[x+2] 
            into the new offsets. offsets[x+1] is not relevant any more */
            new_offsets.push_back(offsets[x]);
            new_offsets.push_back(offsets[x+2]);
            x += 2;
        }
        // if the number of offsets was odd, there might be a dangling offset
        // which we must remember to include in the new_offsets
        if(x+2==offsets.size()) {
            new_offsets.push_back(offsets[x+1]);
        }
        offsets.swap(new_offsets);
    }
    offsets.clear();
    offsets.shrink_to_fit();
}

void SortAndAggregate(std::vector<KmerPairs>& arr) {
    std::vector<KmerPairs>::iterator low,up, it;
    std::vector<KmerPairs> new_arr;

    for( it = arr.begin(); it != arr.end(); ) {
        kmer_t key = (*it).seq;
        low=std::lower_bound(arr.begin(), arr.end(), key, 
            [] (const KmerPairs& lhs, kmer_t rhs) {
                return (lhs.seq < rhs);
        });
 
        up= std::upper_bound (arr.begin(), arr.end(), key,
            [] (kmer_t rhs, const KmerPairs& lhs) {
                return (rhs < lhs.seq);
        });

        
        int sum=0;
        for (auto itr=low; itr!= up; itr++) {
            sum += (*itr).k_count;
        }

        new_arr.push_back(KmerPairs{key, sum});
            it = up;
        }

    arr=new_arr;
    new_arr.clear();
    new_arr.shrink_to_fit();
}

void transfer_kmers (std::vector<int>& scounts_kmer, std::vector<KmerPairs> &kmer_send_buf)  {
    int ssize=0, rsize=0;
    std::vector<int> rcounts_kmer (size,0);
    std::vector<int> rdisp_kmer (size,0);
    std::vector<int> sdisp_kmer (size,0);

    for (int t=0; t<size; t++) ssize += scounts_kmer[t];

    sdisp_kmer[0] = 0;
    for (int i=1; i<size; i++) sdisp_kmer[i] = scounts_kmer[i-1] + sdisp_kmer[i-1];

    //  create contiguous derived data type
    MPI_Datatype rowtype;
    MPI_Type_contiguous(sizeof(KmerPairs), MPI_BYTE, &rowtype);
    MPI_Type_commit(&rowtype);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Alltoall (scounts_kmer.data(), 1, MPI_INT, rcounts_kmer.data(), 1, MPI_INT, MPI_COMM_WORLD);

    for (int t=0; t<size; t++) rsize += rcounts_kmer[t];
    rdisp_kmer[0] = 0;
    for (int i=1; i<size; i++) rdisp_kmer[i] = rcounts_kmer[i-1] + rdisp_kmer[i-1];

    std::vector<KmerPairs> kmer_recv_buf (rsize);
    MPI_Alltoallv(kmer_send_buf.data(), scounts_kmer.data(), sdisp_kmer.data(), rowtype,
                  kmer_recv_buf.data(), rcounts_kmer.data(), rdisp_kmer.data(), rowtype, 
                  MPI_COMM_WORLD);

    kmer_send_buf.clear();
    kmer_send_buf.shrink_to_fit();

    #ifdef FAST 
    ska_sort(kmer_recv_buf.begin(), kmer_recv_buf.end(), [](const auto& i) {return i.seq;});
    #else
    std::sort(kmer_recv_buf.begin(), kmer_recv_buf.end(), [](const auto& i, const auto& j) {return i.seq < j.seq;});
    #endif

    if (kmer_proc_buf.size()) {
        size_t offset = kmer_proc_buf.size();
        kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());

        std::inplace_merge(
            kmer_proc_buf.begin(), 
            kmer_proc_buf.begin() + offset,
            kmer_proc_buf.end(), // this *might* be at the end
                [](const auto& i, const auto& j) { return i.seq < j.seq; }
        );
    } else {
      kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());
    }
    
    kmer_recv_buf.clear();
    kmer_recv_buf.shrink_to_fit();

    SortAndAggregate (kmer_proc_buf);
    
    num_batch_transfers++;
  
    // free datatype
    MPI_Type_free(&rowtype);
    
    /* free the memory */
    rcounts_kmer.clear();
    rdisp_kmer.clear();
    sdisp_kmer.clear(); 
}


void recalculate_min_lmer (kmer_t kmer_in, lmer_t *m_lmer, lmer_t *m_lmer_freq, int *m_pos) {
    lmer_t min_lmer=0, tmp_lmer=0;
    lmer_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, k=0;

    for (k=0; ((KMER_LENGTH-1) - k) >= (LMER_LENGTH-1); k++) {
        lmer_t lmer_out=0;
        for(int j=k; j<LMER_LENGTH+k; j++) {
            lmer_out = kmer_to_lmer (kmer_in, j, lmer_out);
        }

        tmp_lmer = lmer_out;
        tmp_lmer_freq = global_lmer_frequency[tmp_lmer];

        if (k == 0) {
            min_lmer = tmp_lmer;
            min_lmer_freq = tmp_lmer_freq;
            min_pos = 0;
        }
        else {
           if (tmp_lmer_freq < min_lmer_freq) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
               min_pos = k;
           }
        }
    }
    assert (k == (KMER_LENGTH-LMER_LENGTH+1));

    *m_lmer = min_lmer;
    *m_lmer_freq = min_lmer_freq;
    *m_pos = min_pos;
}


void Sliding_window_l (const char *ptr, size_t length) {
  size_t p=0;
  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  lmer_t kmer = 0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+LMER_LENGTH > length) break; /*too short a read*/
    kmer = 0;
    int i;
    for(i=0; ptr[p]!='\n' && i<LMER_LENGTH-1; i++) {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));
    }

    while(p<length && ptr[p]!='\n') {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));
      lmer_frequency[kmer]++;
    }
    p++; /*skip the newline*/
  }
}

void Sliding_window (const char *ptr, size_t length, int *n_kmers, 
    std::vector<std::vector<kmer_t>> &partial_kmer_counts) {

    size_t p=0;
    std::vector<int> scounts_kmer (size,0);
    std::vector<KmerPairs> kmer_send_buf;
    size_t kpos=0;

    /*find start of a read*/
    for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

    int num_kmers=*n_kmers;
    kmer_t kmer = 0; 
    lmer_t lmer_out = 0, min_lmer=0, tmp_lmer=0;
    uint64_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, tmp_pos=0;

    while(p<length) {
        /*skip till newline*/
        for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
        p++; /*skip the newline*/

        if(p+KMER_LENGTH > length) break; /*too short a read*/

        kmer=0, lmer_out=0;
        min_lmer=0, min_lmer_freq=0;
        tmp_lmer=0, tmp_lmer_freq=0;
        min_pos=0, tmp_pos=0;
        int i;

        for(i=0; ptr[p]!='\n' && i<KMER_LENGTH-1; i++) {
            kmer = kmer_shift(kmer, char_to_el(ptr[p]));

            if (i<LMER_LENGTH-1) { 
                lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
            } else {
            lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));

            tmp_lmer = lmer_out;
            tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
            tmp_pos = i-(LMER_LENGTH-1);

            if (i == LMER_LENGTH-1) {
                min_lmer = tmp_lmer;
                min_lmer_freq = tmp_lmer_freq;
            }
            else if (tmp_lmer_freq < min_lmer_freq) {
                min_lmer = tmp_lmer;
                min_lmer_freq = tmp_lmer_freq;
                min_pos = tmp_pos;
            }
        }
        p++;
        }
        
        while(p<length && ptr[p]!='\n') {
            kmer = kmer_shift(kmer, char_to_el(ptr[p]));
            lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
            uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];
    
            if (min_pos < 0) {
                recalculate_min_lmer(kmer, &min_lmer, &min_lmer_freq, &min_pos);
                num_recalculate_lmer++;
            }

            if (lmer_out_freq < min_lmer_freq) {
                min_lmer = lmer_out;
                min_lmer_freq = lmer_out_freq;
                min_pos = KMER_LENGTH-LMER_LENGTH;
            }
            p++;
            min_pos--;

            partial_kmer_counts[retrieve_proc_id(min_lmer)].push_back(kmer);
            num_kmers++;

            if (num_kmers > MAX_KMER_COUNT){
                /* initiate collective communication to pass k-mers and their respective counts to rightful owners
                calculate global owner of each k-mer and populate the k-mer and count to their respective 'p'th local buffers
                if global position of a k-mer in my k_map != me, delete k-mer from my k_map
                iterate through k-mers recieved after collective communication ends, and add k-mers to my k_map 
                reset num_kmers count to 0. */
            
                #ifdef NOOPT
                
                for (int t=0; t<size; t++) {
                    int counter=1;
                    
                    #ifdef FAST
                    ska_sort(kmers_per_proc[t].begin(), kmers_per_proc[t].end());
                    #else 
                    sort(kmers_per_proc[t].begin(), kmers_per_proc[t].end());
                    #endif

                    kmer_t prev=kmers_per_proc[t][0];
                    for(int i = 1; i < (int)(kmers_per_proc[t].size()); i++) {
                        if (kmers_per_proc[t][i] == prev) {
                            counter++;
                        } else {
                            kmer_cnt_tmp_buf[t].push_back(counter);
                            counter=1;
                            prev=kmers_per_proc[t][i];
                        }
                    }
                    
                    kmer_cnt_tmp_buf[t].push_back(counter);

                    kmers_per_proc[t].erase( unique( kmers_per_proc[t].begin(), kmers_per_proc[t].end() ), kmers_per_proc[t].end() );
                    assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
                    scounts_kmer[t] = kmer_cnt_tmp_buf[t].size(); 
                }
                
                #else
                
                for (int t=0; t<size; t++) {
                    int counter=1; kpos=0;

                    #ifdef FAST
                    ska_sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                    #else
                    sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                    #endif

                    kmer_t prev=partial_kmer_counts[t][0];
                    for(int i = 1; i < (int)(partial_kmer_counts[t].size()); i++) {
                        if (partial_kmer_counts[t][i] == prev) {
                            counter++;
                        } else {
                            kmer_send_buf.push_back(KmerPairs{prev, counter});
                            kpos++;
                            counter=1;
                            prev=partial_kmer_counts[t][i];
                        }   
                    }

                    kmer_send_buf.push_back(KmerPairs{prev, counter});
                    kpos++;
                    scounts_kmer[t] = kpos;
                    partial_kmer_counts[t].clear();
                    partial_kmer_counts[t].shrink_to_fit();
                }
            
                #endif
            
                #ifdef NOOPT
                transfer_kmers (scounts_kmer, kmers_per_proc, kmer_cnt_tmp_buf);
                #else
                transfer_kmers (scounts_kmer, kmer_send_buf);
                #endif
                num_kmers = 0;

                #ifdef NOOPT
                for (int t=0; t<size; t++) {
                kmers_per_proc[t].clear();
                kmer_cnt_tmp_buf[t].clear();
                }
                #endif
                
                scounts_kmer.clear(); 
            
            } // end of if condition
        } // end of while loop
        p++; /*skip the newline*/
    }

    *n_kmers = num_kmers;
    scounts_kmer.shrink_to_fit();
}

void process_remaining_kmers(std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
    std::vector<int> scounts_kmer (size,0);
    std::vector<KmerPairs> kmer_send_buf;
    size_t kpos=0;

    for (int t=0; t<size; t++) {
        if(partial_kmer_counts[t].size()) {
            int counter=1;
            kpos=0;

            #ifdef FAST
            ska_sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
            #else
            sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
            #endif

            kmer_t prev=partial_kmer_counts[t][0];
            for(int i = 1; i < (int)(partial_kmer_counts[t].size()); i++) {
                if (partial_kmer_counts[t][i] == prev) {
                    counter++;
                } else {
                    kmer_send_buf.push_back(KmerPairs{prev, counter});
                    kpos++;
                    counter=1;
                    prev=partial_kmer_counts[t][i];
                }		     
            }
            
            kmer_send_buf.push_back(KmerPairs{prev, counter});
            kpos++;
            scounts_kmer[t] = kpos;
            partial_kmer_counts[t].clear();
            partial_kmer_counts[t].shrink_to_fit();
        }
    }

    //clear the partial counts
    for(int k=0; k<size; k++)
        partial_kmer_counts[k].clear();         

    transfer_kmers (scounts_kmer, kmer_send_buf);
    
    scounts_kmer.clear();
    scounts_kmer.shrink_to_fit();
}

void free_kmer_count_buffers() {
    kmer_proc_buf.clear();
    kmer_proc_buf.shrink_to_fit();
}

void perform_kmer_counting (const char *read_data, size_t length) {
    //calculate frequencies of all l-mers in the Read dataset
    Sliding_window_l(read_data, length);

    // Perform Allreduce to obtain global lmer counts across all proc's
    int num_lmers = pow(4, LMER_LENGTH);

    MPI_Allreduce(lmer_frequency.data(), global_lmer_frequency.data(), num_lmers, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    #ifdef NOOPT
    
    std::vector< std::vector<kmer_t> > kmers_per_proc; //, std::vector<kmer_t>);
    std::vector< std::vector<int> > kmer_cnt_tmp_buf; //, std::vector<int>);
    for (int i=0; i<size; i++) {
         kmers_per_proc.push_back(std::vector<kmer_t> ());
         kmer_cnt_tmp_buf.push_back(std::vector<int> ());
    }
    
    #else
    
    std::vector< std::vector<kmer_t> > partial_kmer_counts(size);
    
    #endif
    
    int num_kmers = 0;
   
    #ifdef NOOPT
    Sliding_window(rdata.read_data, rdata.read_data_size, &num_kmers, kmers_per_proc, kmer_cnt_tmp_buf);
    #else
    Sliding_window(read_data, length, &num_kmers, partial_kmer_counts);
    #endif
        
    // initiate communication for the residual kmers
    if (num_kmers) {         
        #ifdef NOOPT
        process_remaining_kmers(kmers_per_proc, kmer_cnt_tmp_buf);
        #else
        process_remaining_kmers(partial_kmer_counts);
        #endif
    }

    #ifdef NOOPT
    
    kmers_per_proc.clear();
    kmers_per_proc.shrink_to_fit();
    kmer_cnt_tmp_buf.clear();
    kmer_cnt_tmp_buf.shrink_to_fit();
    
    #else
    
    for (int i=0; i<size; i++)
        partial_kmer_counts[i].clear();
    partial_kmer_counts.shrink_to_fit(); 
    
    #endif

    lmer_frequency.clear();
    global_lmer_frequency.clear();
}