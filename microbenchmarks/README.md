# Benchmarks to obtain $C_{node}$ and $B_{mem}$ parameters 

## Measuring $C_{node}$

### Compile and Execution
Ensure the processor supports `AVX512` extensions.
```
g++ -std=c++20 -O3 -march=native int64max.cpp -o int64max
./int64max
```
### Sample Output
```
Time taken: 6.0045 seconds
Performance: 121.979 Giga Ops
```
The `Performance` metric tells us the theoretical maximum number of `uint64_t` additions that the current node can perform, and chosen as $C_{node}$ for our analytical model.

## Measuring $B_{mem}$

### Compile and Execution
```
mpicxx -std=c++11 -O3 -o membandwidth membandwidth.cpp
srun -n <total_cores_per_node> membandwidth
```
### Sample Output
```
Read Bandwidth: 47.0343 GB/s
Write Bandwidth: 99.3075 GB/s
```
`Read Bandwidth` is picked as the $B_{mem}$ (memory bandwidth) for our analytical model.