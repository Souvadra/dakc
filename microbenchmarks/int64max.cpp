#include <iostream>
#include <vector>
#include <x86intrin.h>
#include <chrono>
#include <immintrin.h>

void dense_scalar(const uint64_t* vec1, const uint64_t* vec2, uint64_t* out, size_t len) {
  size_t i = 0;
  for (; i + 8 <= len; i += 8) {
    __m512i v1 = _mm512_loadu_si512(&vec1[i]);
    __m512i v2 = _mm512_loadu_si512(&vec2[i]);
    __m512i result = _mm512_add_epi64(v1, v2);
    _mm512_storeu_si512(&out[i], result);
  }
  /* We assume the input length will always be multiple of 8 */
}

int main() {
  const size_t TRIALS = 32000000;       // Number of times to run the test
  const size_t len = 1024;              // Must always be a multiple of 8
  const size_t cores_per_node = 24;     // Number of cores per node (Machine dependent)
  std::vector<uint64_t> vec1(len, 0xFFFFFFFFFFFFFFFF);
  std::vector<uint64_t> vec2(len, 0xFFFFFFFFFFFFFFFF);
  std::vector<uint64_t> out(len);

  volatile size_t trials = TRIALS; // Prevent optimization

  double start = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
  for (size_t i = 0; i < trials; i++) {
    dense_scalar(vec1.data(), vec2.data(), out.data(), len);
  }
  double end = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();

  double time = (end - start) / 1000000000.0;
  std::cout << "Time taken: " << time << " seconds" << std::endl;

  double ops = len * trials;
  double scaling = static_cast<double>(cores_per_node);
  double ops_per_cycle = scaling * (ops / time) / (1024 * 1024 * 1024);
  std::cout << "Max INT64 Performance: " << ops_per_cycle << " Giga Ops" <<  std::endl;

  return 0;
}

// Compiler must support AVX512
// Compile instrution: g++ -std=c++20 -O3 -march=native int64max.cpp
