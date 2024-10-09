# Synthetic Dataset Generation

## How to generate synthetic datasets
Modify the `#define N 5` from `5` to any other variable to generate that many reference synthetic genomes. `N 12` will generate all synthetic references from `Synthetic 20` to `Synthetic 32`.

```
g++ -std=c++17 -O3 inputgenerator.cpp -o step1
./step1
./synthetic_generator.sh
```
All the Synthetic files will be inside the `fqfiles` directory.

## Explanation of each program / script 
- `inputgenerator.cpp` generates `Synthetic 20` to `Synthetic (20+N)` ground truth genomes (`.fa` files). Each genome generated as output of this program is a uniform random string of $2^{20+x}$ length sampled from the $\Sigma = \{A,T,C,G\}$ alphabet $\forall x \in {1, \cdots, N}$.

- `synthetic_generator.sh` uses ART Illumina simulator and generate Illumina like `.fq` files. Then it converts those into into `.fa` format (needed by HySortK and PakMan), and header removed `.txt` format (needed by DAKC).