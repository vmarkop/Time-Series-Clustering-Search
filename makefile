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

CCFLAGS = -lm -g -O2

all: search lsh.out cube.out cluster.out  

search: $(SRCSEARCH)/lsh.cpp
	g++ $(CCFLAGS) -o bin/search $(SRCSEARCH)/lsh.cpp $(SRCLSH)/hashTable.cpp $(SRCSEARCH)/lshUtils.cpp $(SRCLSH)/lsh_frechet_cont.cpp $(SRCLSH)/lsh_frechet_dsc.cpp $(LIB)$(SRCUTIL)/mathUtils.cpp $(LIB)$(SRCUTIL)/projectUtils.cpp $(SRCFRED)/config.cpp $(SRCFRED)/curve.cpp $(SRCFRED)/frechet.cpp $(SRCFRED)/interval.cpp $(SRCFRED)/point.cpp $(SRCFRED)/simplification.cpp -I $(INCSEARCH) -I $(LIB)$(INCUTIL) -I $(INCLSH) -I $(INCFRED)

# lsh.out: $(SRCLSH)/lsh.cpp
# 	g++ $(CCFLAGS) -o bin/lsh.out $(SRCLSH)/lsh.cpp $(SRCLSH)/hashTable.cpp $(SRCLSH)/lshUtils.cpp $(LIB)$(SRCUTIL)/mathUtils.cpp $(LIB)$(SRCUTIL)/projectUtils.cpp -I $(LIB)$(INCUTIL) -I $(INCLSH)

# cube.out: $(SRCCUBE)/cube.cpp
# 	g++ $(CCFLAGS) -o bin/cube.out $(SRCCUBE)/cube.cpp $(SRCCUBE)/HChashTable.cpp $(SRCCUBE)/hypercubeUtils.cpp $(SRC)$(SRCUTIL)/mathUtils.cpp $(SRC)$(SRCUTIL)/projectUtils.cpp -I $(LIB)$(INCUTIL) -I $(INCCUBE)

cluster.out: $(SRCCLUSTER)/cluster.cpp
	g++ $(CCFLAGS) -o bin/cluster.out $(SRCCLUSTER)/cluster.cpp $(SRCSEARCH)/lshUtils.cpp $(SRCLSH)/lsh_frechet_cont.cpp $(SRCLSH)/lsh_frechet_dsc.cpp $(SRCCLUSTER)/clusterUtils.cpp $(SRCLSH)/hashTable.cpp $(SRCCUBE)/HChashTable.cpp $(SRCCLUSTER)/kMeans.cpp $(SRCCLUSTER)/methods.cpp $(SRC)/mathUtils.cpp $(SRC)/projectUtils.cpp -I $(INC) -I $(LIB)$(INCUTIL) -I $(INCCLUSTER) -I $(INCLSH) -I $(INCCUBE) -I $(INCSEARCH)

clean:
	rm -r bin/* outputFiles/*