INC=./lib
INCLSH=./lsh/inc
INCCUBE=./hypercube/inc
INCCLUSTER=./cluster/inc
SRC=./lib
SRCLSH=./lsh/src
SRCCUBE=./hypercube/src
SRCCLUSTER=./cluster/src


CCFLAGS = -lm -g -O3

all: lsh.out cube.out cluster.out  

lsh.out: $(SRCLSH)/lsh.cpp
	g++ $(CCFLAGS) -o bin/lsh.out $(SRCLSH)/lsh.cpp $(SRCLSH)/hashTable.cpp $(SRCLSH)/lshUtils.cpp $(SRC)/mathUtils.cpp $(SRC)/projectUtils.cpp -I $(INC) -I $(INCLSH)

cube.out: $(SRCCUBE)/cube.cpp
	g++ $(CCFLAGS) -o bin/cube.out $(SRCCUBE)/cube.cpp $(SRCCUBE)/HChashTable.cpp $(SRCCUBE)/hypercubeUtils.cpp $(SRC)/mathUtils.cpp $(SRC)/projectUtils.cpp -I $(INC) -I $(INCCUBE)

cluster.out: $(SRCCLUSTER)/cluster.cpp
	g++ $(CCFLAGS) -o bin/cluster.out $(SRCCLUSTER)/cluster.cpp $(SRCCLUSTER)/clusterUtils.cpp $(SRCLSH)/hashTable.cpp $(SRCCUBE)/HChashTable.cpp $(SRCCLUSTER)/kMeans.cpp $(SRCCLUSTER)/methods.cpp $(SRC)/mathUtils.cpp $(SRC)/projectUtils.cpp -I $(INC) -I $(INCCLUSTER) -I $(INCLSH) -I $(INCCUBE)

clean:
	rm -r bin/* outputFiles/*