#ifndef _CLUSTER_UTILS_H_
#define _CLUSTER_UTILS_H_
#include <string>
#include <vector>
#include "mathUtils.h"
#include "projectUtils.h"

#define DEF_VECTOR_HASH_TABLES 3
#define DEF_VECTOR_HASH_FUNCTIONS 4
#define DEF_MAX_NUM_M_CUBE 10
#define DEF_NUM_CUBE_DIM 3
#define DEF_PROBES 2

#define EXIT_FAIL_INPUT_ERR 2
#define EXIT_FAIL_CONFIG_ERR 3
#define EXIT_FAIL_OUTPUT_ERR 4
#define EXIT_FAIL_METHOD_ERR 5

#define CLASSIC_METHOD 0
#define LSH_METHOD 1
#define HYPERCUBE_METHOD 2
#define FRECHET_D_METHOD 3

#define UPDATE_FRECHET 0
#define UPDATE_VECTOR 1

#define TOL (50.0)

typedef struct ClusterStruct *ClusterPtr;
typedef struct ClusterStruct
{
    PointPtr centroidPoint;
    std::vector<PointPtr> points;
    int size;
    double silhouette;
} Cluster;

typedef struct ClusterDataStruct
{

    std::string inputFileName,
        configFileName,
        outputFileName,
        methodName,
        updateName;

    int method, update;

    bool complete;

    int number_of_clusters,
        number_of_vector_hash_tables,
        number_of_vector_hash_functions,
        max_number_M_hypercube,
        number_of_hypercube_dimensions,
        number_of_probes;

    int numberOfInputPoints;
    int dimension;

} clusterInputData;

double silhouette_calculator(PointPtr point,
                             std::vector<Cluster> *clusters,
                             int dimensions);
std::vector<int> get_2_closest_clusters(PointPtr point,
                                        std::vector<Cluster> *clusters,
                                        int dimensions);
PointPtr update_centroid(ClusterPtr cluster, int dimension);
double calculateChanges(std::vector<PointPtr> *centroids, std::vector<Cluster> *clusters,
                        std::vector<PointPtr> *newCentroids,
                        int dimension,
                        int method);
double calculateChangesCurve(std::vector<PointPtr> *centroids,
                             std::vector<Cluster> *clusters,
                             std::vector<PointPtr> *newCentroids,
                             int dimension);
PointPtr computeMeanCurve(PointPtr curve1, PointPtr curve2);
void computeOptimalTraversal(std::vector<std::vector<double>> *_c,
                             std::vector<std::vector<int>> *traversal,
                             int dimension);
PointPtr findMean(treeNodePtr treeNode);
void fillTree(treeNodePtr treeNode,
              std::vector<PointPtr> c_points,
              std::vector<int> *inp);
int getInputData(int argc, char **argv, clusterInputData *CLData);
int execCluster(clusterInputData *CLData,
                std::vector<Cluster> *clusters,
                std::vector<PointPtr> *inputPoints,
                std::vector<PointPtr> *centroidPoints);
double evalSilhouette(clusterInputData *CLData, std::vector<Cluster> *clusters);
int writeToOutput(clusterInputData *CLData,
                  std::vector<Cluster> *clusters,
                  std::vector<PointPtr> *centroidPoints,
                  double totalSilhouette,
                  int tCluster);
void deleteData(std::vector<PointPtr> *inputPoints,
                clusterInputData *CLData,
                std::vector<PointPtr> *centroidPoints);

double minDistBetweenCentroids(std::vector<PointPtr> *centroidPoints, int numOfCentroids, int dimension);

#endif