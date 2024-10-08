# Distributed Asynchronous k-mer counting 

## Input 
Illumina paired-ended or single-ended fastq or fasta files.
For simplicity of reading the fastq files, we preprocess the data to remove the header file, making it simple to read the data in parallel using MPI I/O.

## Directory Structure:
```tree
.
├── dummy_main.cpp (simple code to test the serial parts of the code)
├── src 
│   ├── common 
│   │   └── common.hpp (header file with definitions and classes used by all the kernels)
│   │   └── common.cpp
│   ├── fqreader (read the input fastq/a files, Runtime: MPI + HCLIB Actor)
│   │   ├── fqreader.hpp
│   │   └── fqreader.cpp
│   ├── kcounter (count the k-mers, Runtime: HCLIB Actor)
│   │   ├── kcounter.hpp
│   │   ├── kcounter.cpp
│   └── main
│       ├── parser.hpp (argument parser)
│       └── main.cpp
├── README.md
└── Makefile
```
