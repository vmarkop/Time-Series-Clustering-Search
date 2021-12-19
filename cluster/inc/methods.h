#ifndef _METHODS_H_
#define _METHODS_H_

#include <vector>
#include "../../lsh/inc/hashTable.h"
#include "../../hypercube/inc/HChashTable.h"
#include "clusterUtils.h"
#include "projectUtils.h"
typedef struct duplStruct
{
    PointPtr point;
    std::vector<int> clusters;
} duplicatePoint;

int lloyd_method(std::vector<PointPtr> *centroidPoints, PointPtr point, int dimension);
void lsh_method(HashTables *HashTablesObject, std::vector<PointPtr> *centroids, std::vector<Cluster> *clusters, const std::vector<PointPtr> *inputPoints, clusterInputData *CLData, int numOfInputPoints);
void hyperCube_method(HChashTable *HypercubeObject, std::vector<PointPtr> *centroids, std::vector<Cluster> *clusters, std::vector<PointPtr> *inputPoints, clusterInputData *CLData, int numOfInputPoints);
std::vector<PointPtr> find_duplicates(std::vector<std::vector<PointPtr>> clusterPoints, int numOfClusters);
#endif