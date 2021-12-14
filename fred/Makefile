CC=g++

LSH_OBJS    = main.o point.o simplification.o frechet.o curve.o point.o config.o
CUBE_OBJS = main_cube.o point/point.o util/util.o hashfunctions/hash_func.o hashfunctions/hp.o hashtable_group/Hashtable.o hashtable_group/HashtablesGroup.o hypercube/HypercubeHT.o query_group/Query.o query_group/QueryPointsGroup.o 
CLUSTER_OBJS = main_cluster.o point/point.o util/util.o hashfunctions/hash_func.o hashfunctions/hp.o hashtable_group/Hashtable.o hashtable_group/HashtablesGroup.o hypercube/HypercubeHT.o query_group/Query.o query_group/QueryPointsGroup.o clustering/clustering.o clustering/query_range_set.o

CURVES_LSH_OBJS = main_curve_lsh.o hashtable_curves/HashtableCurves.o  point/point.o util/util.o hashfunctions/hash_func.o hashfunctions/hp.o hashtable_group/Hashtable.o hashtable_group/HashtablesGroup.o hypercube/HypercubeHT.o query_group/Query.o query_group/QueryPointsGroup.o

LFLAGS = -lm -O4 -static-libgcc -static-libstdc++ -ffast-math -std=c++14
METHOD = LSH

INPUT_SMALL = input/input_small_id
QUERY_SMALL = input/query_small_id

OUTPUT_LSH = pro
OUTPUT_CURVES = output/curve_out
OUTPUT_CUBE = output/cube_out
CLUSTER_OUT = output/cluster_out

INPUT_BIG = /tmp/sdi1700262/input_b_id
QUERY_BIG = /tmp/sdi1700262/query_b_id


all: main

curves_lsh: $(CURVES_LSH_OBJS)
	$(CC)  -o $@ $^ $(LFLAGS)


main: $(LSH_OBJS)
	$(CC)  -o $@ $^ $(LFLAGS) 

cube: $(CUBE_OBJS)
	$(CC)  -o $@ $^ $(LFLAGS)

cluster: $(CLUSTER_OBJS)
	$(CC)  -o $@ $^ $(LFLAGS)

%.o: %.c $(HEADER)
	$(CC)  -c -o $@ $< $(LFLAGS)

clean:
	rm -f lsh cube cluster $(LSH_OBJS) $(CUBE_OBJS) $(CLUSTER_OBJS) $(CLUSTER_OUT) $(OUTPUT_LSH) $(OUTPUT_CUBE) $(CURVES_LSH_OBJS) curves_lsh output/curve_out

exec_lsh:
	./lsh -i $(INPUT_SMALL) -q $(QUERY_SMALL) -k 4 -N 1 -R 300 -o $(OUTPUT_LSH) -L 5

exec_curves_lsh:
	./curves_lsh -i $(INPUT_SMALL) -q $(QUERY_SMALL) -k 4 -N 1 -R 300 -o $(OUTPUT_CURVES) -L 5

exec_big_lsh:
	./lsh -i $(INPUT_BIG) -q $(QUERY_SMALL) -k 4 -N 1 -R 150 -o $(OUTPUT_LSH) -L 10


exec_cube:
	./cube -i $(INPUT_SMALL) -q $(QUERY_SMALL) -k 4 -N 1 -R 300 -M 20 -probes 10 -o $(OUTPUT_CUBE) 

exec_cluster:
	./cluster -i $(INPUT_SMALL) -o $(CLUSTER_OUT) -c cluster.conf -m $(METHOD)