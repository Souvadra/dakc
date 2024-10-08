# Actor Sequencing Library - k-mer counting 

Goal of this branch is to match the performance of KMC3 in just counting k-mers. Later, I'll just make it (k+2)-mers and make the cap to 29 (as the max supported k-value) to use it to make de-Bruijn graph out of it.

## Changes to be made in this branch: 
1. We don't need pre_suf branch anymore, we are just counting k-mers initially, until KMC3 performance is matched and beaten [TODO]
2. Instead of stopping after 10k messages, I can make it even more asynchronous, although that'll take some more careful coding. [TODO]

## Input 
Illumina paired-ended or single-ended fastq or fasta files.

## Kernels
Kernel 1: Reach the input files and store them in memory 
Kernel 2: Count the k-mers by parsing through the input data
Kernel 3: Add de-bruijn graph genereation to the k-mer counting
Kernel 4: ...

## Directory Structure:
```tree
.
├── dummy_main.cpp (simple code to test the serial parts of the code)
├── src 
│   ├── common 
│   │   └── common.hpp (header file with definitions and classes used by all the kernels)
│   ├── fqreader (read the input fastq/a files, Runtime: MPI + Selector)
│   │   ├── fqreader.hpp
│   │   └── fqreader.cpp
│   ├── kcounter (count the kmers and generates the de-Bruijn graph, Runtime: Selector)
│   │   ├── kcounter.hpp
│   │   ├── kcounter.cpp
│   │   ├── hash_funcs.c (google's murmurhash3 functions, copied from the HipMer codebase)
│   │   └── hash_funcs.h (same as above)
│   ├── contigs (process the dBG and find the contigs, Runtime: Selector, and a tiny bit MPI (will remove it later))
|   |   ├── contigs.hpp
│   │   ├── contigs.cpp
│   │   └── contighandlers.cpp (the src code for all the message handlers used in this section)
│   ├── falseremover (removes the false positive structures from the de Bruijn graph) [TBD]
│   ├── dummy (set of functions to test different kernels of the code)
│   │   ├── dummy.hpp
│   │   └── dummy.cpp
│   └── main
│       ├── parser.hpp (argument parser)
│       └── main.cpp
├── README.md
└── Makefile
```
