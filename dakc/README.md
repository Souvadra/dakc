# Distributed Asynchronous k-mer counting (DAKC)

## Prerequisite: HCLib Actor Library 

Follow the instructions at `https://hclib-actor.com` to download and install the HCLIB Actor runtime library.

## Input 
Illumina paired-ended or single-ended `FASTQ` files. 
For simplicity of reading the FASTQ files, we preprocess the data to remove the header file, making it simple to read the data in parallel using MPI I/O. 
To avoid confusion, we store these preprocessed files as `.txt`. 
The user can run the `fq2txtmaker.sh` script to generate the header removed input files from a given `FASTQ` file.

### Compile time variables the user should modify based on their use case 
- `KMERLEN`: The length $k$ to use. Current implementation limits $k \leq 32$
- `READLEN`: Length of each read in the input `FASTQ` file.
- `HITTER`: If `HITTER == 0`, then the $L_3$ aggregation protocol is not performed, and vice versa.
- `BIGKSIZE`: `2 x BIGKSIZE` is the $C_2$ parameter size, mentioned in the paper.
- `KCOUNT_BUCKET_SIZE`: The value of this parameter determines $C_3$ parameter value. 
- `BENCHMARK`: If present, the program will generate statistics regarding the program behavior and output.

## How to compile

Open the Makefile and update `COMPILETIMEVARS` accordingly. 
[Note: The user does not need to change anything for executing `DAKC` on synthetic datasets generated using data generation scripts provided to the user.]

Type `make clean && make` to compile DAKC.

## How to execute 
```
srun -N <num_nodes> -n <total_cores> --cpu-bind=cores dakc -f <input_file>
```

Note: we recommend creating one process per physical core of the CPU for optimal performance. 
In the above `srun` command, `<total_cores>` should be the total number of physical cores present in all the nodes being used for the execution.

## Directory structure:
```tree
.
├── src 
│   ├── common 
│   │   └── common.hpp (header file with definitions and classes used by all the kernels)
│   │   └── common.cpp
│   ├── fqreader (read the input fastq/a files, Runtime: MPI + HCLIB Actor)
│   │   ├── fqreader.hpp
│   │   └── fqreader.cpp
│   ├── kcounter (count the k-mers, Runtime: HCLIB Actor)
│   │   ├── ska_sort.hpp
│   │   ├── kcounter.hpp
│   │   ├── kcounter.cpp
│   └── main
│       ├── parser.hpp (argument parser)
│       └── main.cpp
├── README.md
└── Makefile
```
