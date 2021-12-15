#ifndef _LSH_UTILS_H_
#define _LSH_UTILS_H_

#include <string>
#include "mathUtils.h"
#include "projectUtils.h"
// #include "hashTable.h"

#define DEF_K 4
#define DEF_L 5

#define ALG_LSH 0
#define ALG_HC 1
#define ALG_FR 2

#define MTR_DISC 0
#define MTR_CONT 1

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
} inputData;

inputData *getInputData(int *argc, char **argv);
int writeToOutput(inputData *LSHData,
                  std::vector<PointPtr> queryPoints,
                  std::vector<kNeighboursPtr> queryOutputData,
                  std::vector<kNeighboursPtr> queryTrueNeighbors,
                  double tLSH,
                  double tTrue,
                  std::string algorithm);

int writeToOutputFrDsc(inputData *SearchData,
                       std::vector<PointPtr> queryPoints,
                       std::vector<NeighbourPtr> queryOutputData,
                       std::vector<NeighbourPtr> queryTrueNeighbors,
                       double tLSH,
                       double tTrue,
                       std::string algorithm);

void deleteData(std::vector<PointPtr> *inputPoints,
                std::vector<PointPtr> *queryPoints,
                std::vector<std::vector<Neighbour> *> *k_nearest_neighbours,
                std::vector<kNeighboursPtr> *queryOutputData,
                std::vector<kNeighboursPtr> *queryTrueNeighbors,
                inputData *LSHData);
#endif