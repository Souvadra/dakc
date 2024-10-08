#include <iostream>
#include <sys/time.h>
#include <vector>
#include <unordered_map> // should replace with some other hash map
#include <cstring>
#include <shmem.h>

#ifdef PROFILE
#define ENABLE_TCOMM_PROFILING
#endif

#include "selector.h"

#include <mpi.h> // Just for the input file reading part we'll need MPI 

#include "parser.hpp"
#include "common.hpp"
#include "fqreader.hpp"
#include "kcounter.hpp"

int main(int argc, char** argv) {
    // initialize the MPI runtime 
    int rank, size; // rank and size will be used for MPI_PEs only
    MPI_Init(0,0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // initialize shmem + selector runtime 
    shmem_init();
    const char *deps[] = {"system", "bale_actor"};

    hclib::launch(deps, 2, [=] {
        // parse the command line arguments 
        arg_parser arg(argc, argv);
        if (rank == 0) arg.print_params();

        shmem_barrier_all();
        std::unordered_map<kmer_t, count_t> dbg; // Is there any use for this anymore ??
        std::vector<kmer_t> vectordbg;
        
        // read the fasta/q files (kernel 1, part 1)
        fqreader fq(arg.file_name, READLEN, rank, size);
        char* read_chunk = fq.read_file();
        
        // time to perform k-mer counting 
        kmercounter km(read_chunk, vectordbg);
        
        // free the variables
        free(read_chunk);
    });

    // finalize shmem
    shmem_finalize();
    return 0;
}
