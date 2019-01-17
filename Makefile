project := hyperPraw

CXX = mpicxx
ARCHER_CC = CC
SRC_PATH = ./src
SRC_EXT = cpp

FILES = $(shell find $(SRC_PATH) -name '*.$(SRC_EXT)' | sort -k 1nr |cut -f2-) 
OUT = $(project)

INCLUDES = -I./include -I$(HOME)/Zoltan_v3.83/build/include 
LIBS = -L$(HOME)/Zoltan_v3.83/build/lib
 
INCLUDES_INTEL = -I./include -I$(HOME)/Zoltan_v3.83_intel/build/include 
INTEL_LIBS = -L$(HOME)/Zoltan_v3.83_intel/build/lib 

INCLUDES_ARCHER = -I./include -I$(HOME)/Zoltan_v3.83/build/include 
ARCHER_LIBS = -L$(HOME)/Zoltan_v3.83/build/lib 

LDFLAGS = -lzoltan
LSTATIC = 

build: $(FILES)
	$(CXX) -O3 -o $(OUT) $(INCLUDES) $(LIBS) $(FILES) $(LDFLAGS) $(LSTATIC) --std=c++11

.PHONY: debug
debug:
	$(CXX) -g -o $(OUT) $(INCLUDES) $(LIBS) $(FILES) $(LDFLAGS) $(LSTATIC) --std=c++11

.PHONY: archer
archer:
	$(ARCHER_CC) -O3 -o $(OUT) $(INCLUDES_ARCHER) $(ARCHER_LIBS) $(FILES) $(LDFLAGS) $(LSTATIC)


.PHONY: clean
clean:
	rm $(OUT)

