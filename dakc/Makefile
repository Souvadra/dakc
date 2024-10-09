include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.pre.mak
include $(HCLIB_ROOT)/include/hclib.mak
include $(HCLIB_ROOT)/../modules/bale_actor/inc/hclib_bale_actor.post.mak

CXX ?= CC
SRUN ?= SRUN
BALE_FLAGS ?= -DUSE_SHMEM=1

CFLAGS = -std=c++17 -O3 -march=native $(BALE_FLAGS) $(HCLIB_CFLAGS)
LIBS = $(HCLIB_LDFLAGS) $(HCLIB_LDLIBS) -lspmat -lconvey -lexstack -llibgetput -lhclib_bale_actor -lm -loshmem -lmpi 
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