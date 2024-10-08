include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

CXX ?= CC
SRUN ?= SRUN
BALE_FLAGS ?= -DUSE_SHMEM=1

CFLAGS = -std=c++17 -O3 -march=native $(BALE_FLAGS) $(HCLIB_CFLAGS) # Add the AVX512 flgas here for Intel machines 
LIBS = $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm -loshmem -lmpi 
# -lsohmem and lmpi are not required in Phoenix with GCC compilers
COMPILETIMEVARS = -DKMERLEN=31 -DREADLEN=150 -DMIN_KMER_COUNT=0 -DHITTER=0 -DBIGKSIZE=16 -DKCOUNT_BUCKET_SIZE=10000 # -DBENCHMARK

COMMON = -I$(PWD)/src/common
FQREADER = -I$(PWD)/src/fqreader
KCOUNTER = -I$(PWD)/src/kcounter
MAIN = -I$(PWD)/src/main

INCLUDE = $(COMMON) $(FQREADER) $(KCOUNTER) $(MAIN)

# Compile time variables needed for profiling
PAPI_ROOT=/usr/local/pace-apps/manual/packages/papi/7.0.1/usr/local
PROFILE_FLAGS = -I${PAPI_ROOT}/include -L${PAPI_ROOT}/lib -lpapi -DPROFILE

# Assuming your project is C++, adjust as necessary for C projects
SRC = $(wildcard src/*/*.cpp) 
OBJ = $(SRC:.cpp=.o) 
HDR = $(wildcard src/*/*.hpp) 
PROFILE_OBJ := $(SRC:.cpp=.profile.o)
NAME = dakc

all: $(NAME)

$(NAME): $(OBJ)
	$(CXX) $(CFLAGS) $(INCLUDE) $(COMPILETIMEVARS) -o $(NAME) $^ $(LDFLAGS) $(LIBS)

%.o: %.cpp $(HDR)
	$(CXX) $(CFLAGS) $(INCLUDE) $(COMPILETIMEVARS) -c -o $@ $<

%.profile.o: %.cpp $(HDR)
	$(CXX) $(CFLAGS) $(INCLUDE) $(COMPILETIMEVARS) $(PROFILE_FLAGS) -c -o $@ $<

profile: $(PROFILE_OBJ)
	$(CXX) $(CFLAGS) $(INCLUDE) $(COMPILETIMEVARS) $(PROFILE_FLAGS) -o $(NAME) $^ $(LDFLAGS) $(LIBS)

clean: 
	rm -f $(wildcard src/*/*.o) $(NAME) 

# Scripts to run the program on various datasets  

launch:
	srun -n 16 --cpu-bind=cores dakc -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_1048576/synthetic_1048576.txt

llnode:
	srun -n 24 --cpu-bind=cores dakc -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_134217728/synthetic_134217728.txt

ll:
	srun -n 192 --cpu-bind=cores dakc -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_134217728/synthetic_134217728.txt

NN ?= 24
# Run Pakman on synthetic data
S20:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_1048576/synthetic_1048576.txt

S21:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_2097152/synthetic_2097152.txt

S22:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_4194304/synthetic_4194304.txt

S23:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_8388608/synthetic_8388608.txt

S24:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_16777216/synthetic_16777216.txt

S25: 
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_33554432/synthetic_33554432.txt

S26:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_67108864/synthetic_67108864.txt

S27:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_134217728/synthetic_134217728.txt

S28:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_268435456/synthetic_268435456.txt

S29:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_536870912/synthetic_536870912.txt

S30:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_1073741824/synthetic_1073741824.txt

S31:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_2147483648/synthetic_2147483648.txt

S32:
	srun -n $(NN) --cpu-bind=cores wf3 -f /storage/coda1/p-vsarkar9/1/shati3/synthetic_genomes/fqfiles/synthetic_4294967296/synthetic_4294967296.txt

# Run on real data
R1: # Phocoena sinus
	echo "Make srue you recompiled with READLEN 151"
	srun -n $(NN) --cpu_bind=cores wf3 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR25743144/SRR25743144_1.txt

R2: # F. vesca
	srun -n $(NN) --cpu_bind=cores wf3 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR26113965/SRR26113965_1.txt

R3: # H. sapiens
	echo "Make srue you recompiled with READLEN 149"
	srun -n $(NN) --cpu_bind=cores wf3 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR28206931/SRR28206931_1.txt

R4: # S. coelicolor
	srun -n $(NN) --cpu_bind=cores wf3 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR28892189/SRR28892189_1.txt

R5: # P. aeruginosa
	srun -n $(NN) --cpu_bind=cores wf3 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR29163078/SRR29163078_1.txt

R6: # Wheat 
	srun -n $(NN) --cpu_bind=cores wf3 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR29871703/SRR29871703_1.txt

R7: # Salamander 
	echo "Make srue you recompiled with READLEN 125"
	srun -n $(NN) --cpu_bind=cores wf3 -f /storage/home/hcoda1/6/shati3/scratch/real_genomics/SRR7443702/SRR7443702_1.txt
