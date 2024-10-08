#include <iostream> 
#include <string> 
#include <sstream>
#include <vector> 
#include <cstring>
#include <cmath>
#include <cassert>

#include <mpi.h>

#include "fqreader.hpp"
#include "common.hpp"

char* fqreader::read_file() {
    uint64_t numreads = 0, total_reads = 0, counted_reads = 0;
    size_t read_data_size = 0;
    double starttime, endtime, local_readingtime, global_readingtime;
    char *chunk; 
    MPI_Offset filesize, start, end, readbuf;

    if (!is_txt) {
      std::cout << "FA and FQ files are not natively supported !!" << std::endl;
    }

    if (rank == 0) 
      std::cout << "Start reading the input dataset(s)" << std::endl;

    int ierr = MPI_File_open(MPI_COMM_WORLD, filename.c_str(), 
      MPI_MODE_RDONLY, MPI_INFO_NULL, &inputfile);
    
    if (ierr) {
      if (rank == 0) 
          std::cout << "Could not open the input (modified) FASTA/Q file" << std::endl;
      MPI_Finalize();
      exit(2);
    }

    // Start the timer
    starttime = MPI_Wtime();

    MPI_File_get_size(inputfile, &filesize);
    
    total_reads = (filesize + 1) / (read_len + 1); // each character is 1 BYTE
    localsize = total_reads / size;
    localsize *= (read_len + 1);

    start = rank * localsize;
    end = start + localsize - 1;

    if (rank == size - 1) end = (filesize - 1);

    localsize = end - start + 1;

  #if DEBUG 
    std::cout << "PE: " << start << ", " << end << std::endl;
  #endif

    // Provide enough space for the string storing the local data
    chunk = (char*) malloc( (localsize + 1) * sizeof(char) );
    chunk[localsize] = '\0';

    // Set the file view for each process
    MPI_File_set_view(inputfile, start, MPI_CHAR, MPI_CHAR, "native", 
                     MPI_INFO_NULL);

    // Read the chunk of the file
    MPI_File_read(inputfile, chunk, localsize, MPI_CHAR, MPI_STATUS_IGNORE);
    // MPI_Barrier(MPI_COMM_WORLD); // should I used MPI_File_read_all in this step? 
    endtime = MPI_Wtime();

    // free the variables 
    MPI_File_close(&inputfile);

    // calculate the time taken in reading the inputs
    local_readingtime = endtime - starttime;
    MPI_Barrier(MPI_COMM_WORLD); // Is this even required ?
    MPI_Reduce(&local_readingtime, &global_readingtime, 1, MPI_DOUBLE, 
               MPI_MAX, 0, MPI_COMM_WORLD);

    MPI_Reduce(&numreads, &counted_reads, 1, MPI_UINT64_T, MPI_SUM,
               0, MPI_COMM_WORLD);

    // Print the performance metric
    if (rank == 0) {
        std::cout << "Reading time: " << global_readingtime
        << " seconds using " << size << " processors." << std::endl;
    }

    return chunk;
}