# Asynchronous Distributed Memory Parallel $k$-mer counting Algorithm

This is the global repository of all algorithms and software developed for the IPDPS 2025 article on asynchronous, distributed-memory $k$-mer counting algorithm. 
Each directory contains a separate README file describing how to compile, modify, and execute each software. 

**Link to the article**: <will be provided once IPDPS releases its 2025 conference proceedings> 

## Directory Structure

- **/dakc**: Distributed-memory asynchronous $k$-mer counting implementation of `DAKC` algorithm.
- **/analytical_model**: Analytical model for distributed-mmemory $k$-mer counting. 
- **/microbenchmarks**: Benchmark programs used to obtain system parameters for modeling asynchronous $k$-mer counting on supercomputers.
- **/pakman_star**: Faster version of PakMan algorithm (obtained from the original GitHub repository), a state-of-the-art distributed memory short-read genome assembly toolkit.
- **/syngen**: Software to generate synthetic datasets used in the article. 
