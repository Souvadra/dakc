#ifndef FQREADER_H
#define FQREADER_H

#include <iostream> 
#include <string> 
#include <vector>
#include <sstream>

#include <mpi.h>

#include "common.hpp"

class fqreader { 
public: 
    MPI_File inputfile;
    bool is_fq = true;
    bool is_txt = true; 
    int rank, size, read_len; 
    std::string filename;
    MPI_Offset localsize;

    fqreader(std::string filename, const int read_length, const int rank, 
            const int size) { 
        this->rank = rank; 
        this->size = size; 
        this->read_len = read_length;
        this->filename = filename;

        // break the filename based on delimeter '.'
        std::istringstream iss(this->filename);
        std::string token;

        while (std::getline(iss, token, '.')) { 
            // std::cout << token << std::endl;
            // do nothing
        } 

        // update the is_fq flag depending upon the token
        if (token != "fq") this->is_fq = false;
        if (token != "txt") this->is_txt = false;

        // opportunity to serially peek into the file and 
        // get the readbuf information dynamically
    }

    char* read_file();

private: 
    bool saw_at = false;
    bool saw_plus = false;
};

#endif