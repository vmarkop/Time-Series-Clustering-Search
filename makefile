INCLSH=./lsh/inc
INCCUBE=./hypercube/inc
INCCLUSTER=./cluster/inc
INCUTIL=/utils/inc
LIB=./lib
SRCLSH=./lsh/src
SRCCUBE=./hypercube/src
SRCCLUSTER=./cluster/src
SRCUTIL=/utils/src

SRCSEARCH=./search/src
INCSEARCH=./search/inc

SRCFRED=./lib/Fred/src
INCFRED=./lib/Fred/include

INCPY=./lib/pybind

SRCUNIT=./unit

CCFLAGS = -lm -g -O3
UNITFLAGS = -lm -lpthread -lgtest

all: search cluster unit  

search: $(SRCSEARCH)/search.cpp
	g++ $(CCFLAGS) -o bin/search.out $(SRCSEARCH)/search.cpp $(SRCLSH)/hashTable.cpp $(SRCCUBE)/HChashTable.cpp $(SRCSEARCH)/lshUtils.cpp $(SRCLSH)/lsh_frechet_cont.cpp $(SRCLSH)/lsh_frechet_dsc.cpp $(LIB)$(SRCUTIL)/mathUtils.cpp $(LIB)$(SRCUTIL)/projectUtils.cpp $(SRCFRED)/config.cpp $(SRCFRED)/curve.cpp $(SRCFRED)/frechet.cpp $(SRCFRED)/interval.cpp $(SRCFRED)/point.cpp $(SRCFRED)/simplification.cpp -I $(INCSEARCH) -I $(LIB)$(INCUTIL) -I $(INCLSH) -I $(INCFRED) -I $(INCCUBE)

cluster: $(SRCCLUSTER)/cluster.cpp
	g++ $(CCFLAGS) -o bin/cluster.out $(SRCCLUSTER)/cluster.cpp $(SRCSEARCH)/lshUtils.cpp $(SRCLSH)/lsh_frechet_dsc.cpp $(SRCCLUSTER)/clusterUtils.cpp $(SRCLSH)/hashTable.cpp $(SRCCUBE)/HChashTable.cpp $(SRCCLUSTER)/kMeans.cpp $(SRCCLUSTER)/methods.cpp $(LIB)$(SRCUTIL)/mathUtils.cpp $(LIB)$(SRCUTIL)/projectUtils.cpp -I $(LIB)$(INCUTIL) -I $(INCCLUSTER) -I $(INCLSH) -I $(INCCUBE) -I $(INCSEARCH)

unittest: $(SRCUNIT)/unit_testing.cpp
	g++ -o bin/unit_test.out $(SRCUNIT)/unit_testing.cpp $(SRCLSH)/hashTable.cpp $(SRCCUBE)/HChashTable.cpp $(SRCSEARCH)/lshUtils.cpp $(SRCLSH)/lsh_frechet_cont.cpp $(SRCLSH)/lsh_frechet_dsc.cpp $(LIB)$(SRCUTIL)/mathUtils.cpp $(LIB)$(SRCUTIL)/projectUtils.cpp $(SRCCLUSTER)/clusterUtils.cpp $(SRCCLUSTER)/kMeans.cpp $(SRCCLUSTER)/methods.cpp $(SRCFRED)/* -I $(INCSEARCH) -I $(LIB)$(INCUTIL) -I $(INCLSH) -I $(INCFRED) -I $(INCCUBE) -I $(INCCLUSTER) $(UNITFLAGS)

clean:
	rm -r bin/* outputFiles/*