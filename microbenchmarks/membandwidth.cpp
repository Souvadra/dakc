#include <iostream>
#include <chrono>
#include <vector>
#include <utility>

#include <mpi.h>

const size_t ARRAY_SIZE = 1000000000;  // Size of the array (100 million elements)

std::pair<double, double> measure_bandwidth() {
  std::vector<int> array(ARRAY_SIZE, 0);

  // Write to the array to ensure it is allocated in DRAM
  for (size_t i = 0; i < ARRAY_SIZE; ++i) {
    array[i] = i;
  }

  // Measure read bandwidth
  auto start = std::chrono::high_resolution_clock::now();
  volatile int sum = 0;
  for (size_t i = 0; i < ARRAY_SIZE; ++i) {
    sum += array[i];
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = end - start;
  double read_bandwidth = (ARRAY_SIZE * sizeof(int)) / (duration.count() * 1e9);  // GB/s

  // Measure write bandwidth
  start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < ARRAY_SIZE; ++i) {
    array[i] = i;
  }
  end = std::chrono::high_resolution_clock::now();
  duration = end - start;
  double write_bandwidth = (ARRAY_SIZE * sizeof(int)) / (duration.count() * 1e9);  // GB/s

  return std::pair<double, double>(read_bandwidth, write_bandwidth);
}

int main() {
  /* initialize MPI */
  MPI_Init(NULL, NULL);
  
  int rank, size;
  double lreadbw, lwritebw, greadbw, gwritebw;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::pair<double, double> bwpair = measure_bandwidth();

  lreadbw = bwpair.first;
  lwritebw = bwpair.second;

  MPI_Reduce(&lreadbw, &greadbw, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&lwritebw, &gwritebw, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "Read Bandwidth: " << greadbw << " GB/s" << std::endl;
    std::cout << "Write Bandwidth: " << gwritebw << " GB/s" << std::endl;
  }

  MPI_Finalize();

  return 0;
}

// Compiliation: mpicxx -std=c++11 -O3 -o membandwidth membandwidth.cpp