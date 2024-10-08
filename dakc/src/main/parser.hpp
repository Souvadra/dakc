#ifndef PARSER_HPP
#define PARSER_HPP

#include <iostream>
#include <string>
#include <assert.h>

#include <mpi.h>
#include <getopt.h> // for argument parsing

#include "common.hpp"

// table of all supported options in their long form 
option longopts[] { 
  {"help", no_argument, NULL, 'h'},
  {"file", required_argument, NULL, 'f'}, 
  {0}
};

class arg_parser { 
public:
  // arguments of the program
  std::string     file_name = "0";

  // description of al supported options
  void print_usage();

  // make sure the argument parsing has performed succesfully
  void arg_parser_sanity_check();

  // print all the parameters once reading is done
  void print_params();

  // default and only constructor                                            
  arg_parser(const int argc, char** const argv) { 
    bool help_flag = false;
    int opt;

    while((opt = getopt_long(argc, argv, "hp:f:g:r:k:b:m:x:z:y:", longopts, 0)) != -1) { 
      
      switch (opt) { 
        case 'h':
          print_usage();
          break;
        case 'f':
          this->file_name.assign(optarg);
          break;
        default:
          print_usage();
          assert(0 && "Should not reach here !!");
      }
    }

    arg_parser_sanity_check();
  }
};

inline void arg_parser::print_usage() { 
  std::cout << "required program parameters:" << std::endl;
  std::cout << "-h, --help\t" << "Print this help and exit" << std::endl;
  std::cout << "-f, --file1\t" << "file name" << std::endl;
}

inline void arg_parser::arg_parser_sanity_check() { 
  // Must provide file name
  assert(this->file_name != "0");
}

inline void arg_parser::print_params() {
  std::cout << "File Name : " << this->file_name << std::endl; 
  std::cout << "Read Length : " << READLEN << std::endl;
  std::cout << "k-mer Length : " << KMERLEN << std::endl;
  std::cout << "C3 Length : " << KCOUNT_BUCKET_SIZE << std::endl;
  std::cout << "C2 Length : " << BIGKSIZE * 2 << std::endl;
  // std::cout << "minimizer_length = " << MINIMIZERLEN << std::endl;
  // std::cout << "maximum kmer count = " << MAX_KMER_COUNT << std::endl; 
  // std::cout << "min kmer count = " << MIN_KMER_COUNT << std::endl;
  // std::cout << "min contig len = " << MINCONTIGLEN << std::endl;
}

#endif
