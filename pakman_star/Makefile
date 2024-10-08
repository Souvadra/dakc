CC = mpic++
CFLAGS = -std=c++14 -O3 -fopenmp -DLMER_LENGTH=8  -DFAST # -DPERFORM_CONTIG_GENERATION  -DCOMPACT_PGRAGH
INCLUDE = -I$(PWD)/includes

GET_WINSIZE=0
KVAL=$(shell echo $(ksize)-$(GET_WINSIZE) | bc)

ifeq ($(shell expr $(KVAL) \<= 0), 1)
  KVAL=30
endif

ifeq ($(shell echo $(ksize) \> 64), 1)
$(error ksize is greater than 64)
endif

ifeq ($(shell expr $(KVAL) \<= 31), 1)
  CFLAGS += -DWINDW_SIZE=$(KVAL)
else
  CFLAGS += -DWINDW_SIZE=$(KVAL) -DEXTEND_KMER
endif 

ENABLE_DEBUG=0
ifeq ($(ENABLE_DEBUG), 1)
  LDFLAGS=-lm -ldl -fsanitize=address -g 
else
  LDFLAGS=-lm -ldl
endif
 
ENABLE_DUMP_DEBUG_DATA=0
ifeq ($(ENABLE_DUMP_DEBUG_DATA), 1)
  CFLAGS += -DDEBUG_KSIZE -DDEBUG_WIRE -DDEBUG_IDSET
endif    

SRC = $(wildcard src/*.cpp)
HDR = $(wildcard includes/*.h)
NAME = pakman

all: $(NAME)

$(NAME) : $(SRC) $(HDR)
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) -o $(NAME) $? -fopenmp

clean:
	rm -f $(NAME)

distclean: clean

NN ?= 128
# Run Pakman on synthetic data
S20:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_1048576/synthetic_1048576.fa

S21:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_2097152/synthetic_2097152.fa

S22:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_4194304/synthetic_4194304.fa

S23:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_8388608/synthetic_8388608.fa

S24:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_16777216/synthetic_16777216.fa

S25: 
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_33554432/synthetic_33554432.fa

S26:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_67108864/synthetic_67108864.fa

S27:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_134217728/synthetic_134217728.fa

S28:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_268435456/synthetic_268435456.fa

S29:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_536870912/synthetic_536870912.fa

S30:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_1073741824/synthetic_1073741824.fa

S31:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_2147483648/synthetic_2147483648.fa

S32:
	srun -n $(NN) --cpu-bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_4294967296/synthetic_4294967296.fa

# Run on real data
R1: # Phocoena sinus
	srun -n $(NN) --cpu_bind=cores pakman -r 151 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR25743144/SRR25743144_1.fa .

R2: # F. vesca
	srun -n $(NN) --cpu_bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR26113965/SRR26113965_1.fa .

R3: # H. sapiens
	srun -n $(NN) --cpu_bind=cores pakman -r 149 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR28206931/SRR28206931_1.fa .

R4: # S. coelicolor
	srun -n $(NN) --cpu_bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR28892189/SRR28892189_1.fa .

R5: # P. aeruginosa
	srun -n $(NN) --cpu_bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR29163078/SRR29163078_1.fa .

R6: # Wheat 
	srun -n $(NN) --cpu_bind=cores pakman -r 150 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR29871703/SRR29871703_1.fa .

R7: # Salamander 
	srun -n $(NN) --cpu_bind=cores pakman -r 125 -c 50 -b 1000000000 -t 21 -n 100000 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR7443702/SRR7443702_1.fa .
