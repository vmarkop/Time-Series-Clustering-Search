#ifndef _LSH_UTILS_H_
#define _LSH_UTILS_H_

#include <string>
#include "../../lib/mathUtils.h"
#include "../../lib/projectUtils.h"
// #include "hashTable.h"

#define DEF_K 4
#define DEF_L 5
#define ALG_LSH 0
#define ALG_HC 1
#define ALG_FR 2

typedef struct searchInputDataStruct
{
    std::string inputFileName;
    std::string queryFileName;
    std::string outputFileName;
    int numberOfHyperplanes; // intK numberOfHyperplanes
    int intL;
    int intM;
    int probes;
    int algorithm;
    int metric;
    double delta;
    int dimension;
    bool distance_true_visible;
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
                inputData *LSHData);
#endif