#include <iostream>
#include <vector>
#include <x86intrin.h> // Header for _rdtsc()
#include <chrono>
#include <immintrin.h>

void dense_scalar(const uint64_t* vec1, const uint64_t* vec2, uint64_t* out, size_t len) {
  // Example implementation (replace with actual logic)
  size_t i = 0;
  for (; i + 8 <= len; i += 8) {
    __m512i v1 = _mm512_loadu_si512(&vec1[i]);
    __m512i v2 = _mm512_loadu_si512(&vec2[i]);
    __m512i result = _mm512_add_epi64(v1, v2);
    _mm512_storeu_si512(&out[i], result);
  }
  /* Let's assume the input length will always be multiple of 8 */
  // // Handle remaining elements
  // for (; i < len; i++) {
  //   out[i] = vec1[i] + vec2[i];
  // }
}

int main() {
  const size_t TRIALS = 32000000;
  const size_t len = 1024;
  std::vector<uint64_t> vec1(len, 0xFFFFFFFFFFFFFFFF);
  std::vector<uint64_t> vec2(len, 0xFFFFFFFFFFFFFFFF);
  std::vector<uint64_t> out(len);

  volatile size_t trials = TRIALS; // Prevent optimization

  // uint64_t start1 = _rdtsc();
  double start = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
  for (size_t i = 0; i < trials; i++) {
    dense_scalar(vec1.data(), vec2.data(), out.data(), len);
  }
  double end = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
  // uint64_t end1 = _rdtsc();

  double time = (end - start) / 1000000000.0;
  std::cout << "Time taken: " << time << " seconds" << std::endl;

  // std::cout << "Time taken: " << (end1 - start1) << " time" << std::endl;

  // std::cout << "Time taken: " << (end1 - start1) << " cycles" << std::endl;

  /* I want number of operations done per second printed out */
  // double cycles = end1 - start1;
  double ops = len * trials;
  double ops_per_cycle = (ops / time) / (1024 * 1024 * 1024);
  std::cout << "Performance: " << ops_per_cycle << " Giga Ops" <<  std::endl;

  return 0;
}