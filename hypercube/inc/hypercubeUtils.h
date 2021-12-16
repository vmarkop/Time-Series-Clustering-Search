#ifndef _HYPERCUBE_UTILS_H_
#define _HYPERCUBE_UTILS_H_

#include <string>
#include "HChashTable.h"
#include "mathUtils.h"

#define DEF_K 14
#define DEF_M 10
#define DEF_PROBES 2
#define DEF_N 1
#define DEF_R 500

typedef struct HCinputDataStruct
{
    std::string inputFileName;
    std::string queryFileName;
    std::string outputFileName;
    int projectionDimension;
    int maxCandidatePoints;
    int probes;
    int numberOfNearest;
    int radius;
    int dimension;
} inputData;

inputData *getInputData(int *argc, char **argv);
int writeToOutput(inputData *LSHData,
                  std::vector<PointPtr> queryPoints,
                  std::vector<kNeighboursPtr> queryOutputData,
                  std::vector<kNeighboursPtr> queryTrueNeighbors,
                  std::vector<std::vector<PointPtr>> queryRangeSearch,
                  std::vector<double> tLSH,
                  std::vector<double> tTrue);
void deleteData(std::vector<PointPtr> *inputPoints,
                std::vector<PointPtr> *queryPoints,
                std::vector<std::vector<Neighbour> *> *k_nearest_neighbours,
                std::vector<kNeighboursPtr> *queryOutputData,
                std::vector<kNeighboursPtr> *queryTrueNeighbors,
                inputData *HCData);

#endif