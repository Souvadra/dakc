CC = mpic++
CFLAGS = -std=c++14 -O3 -fopenmp -DLMER_LENGTH=8 -DFAST # -DPERFORM_CONTIG_GENERATION  -DCOMPACT_PGRAGH
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