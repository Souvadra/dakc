#include <iostream> 
#include <random> 
#include <string> 

#define N 2 // Number of files to generate
#define START 20
static const char bases[] = {'A', 'C', 'G', 'T', 'N'};

inline char int2base(int x) {
  return bases[x];
}

int main() {
  // Seed the random number generator
  std::random_device rd;
  std::mt19937 gen(rd());
  
  // Define the distribution for integers between 0 and 3 (inclusive)
  std::uniform_int_distribution<> dis(0, 3);


  uint64_t length = 0, curr_length = 0;
  std::string line, filename;

  for (int i = 0; i < N; i++) {
    
    length = static_cast<uint64_t>(1) << (START + i);
    curr_length = 0;
    filename = "synthetic_" + std::to_string(START+i) + ".fa";
    
    // open a file for writing using the filename 
    FILE *file = fopen(filename.c_str(), "w");  
    if (file == NULL) {
      std::cerr << "Error: unable to open file " << filename << std::endl;
      return 1;
    }

    // Write the header
    fprintf(file, ">synthetic_%llu\n", length);

    // Start writing the sequence 
    curr_length = 0;

    while (curr_length < length) {
      line = "";
      line.reserve(70);

      for (int k = 0; k < 70; k++) {
        line += int2base(dis(gen));
        curr_length++;
        if (curr_length == length) {
          fprintf(file, "%s\n", line.c_str());
          break;
        }
      }
      if (curr_length == length) {
        break;
      }
      fprintf(file, "%s\n", line.c_str());
    }

  }
  return 0;
}

// Compile instruction: g++ -std=c++17 -O3 inputgenerator.cpp
