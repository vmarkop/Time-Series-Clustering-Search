#ifndef _HASH_TABLE_H_
#define _HASH_TABLE_H_

#include <stdlib.h>
#include <vector>
#include <string>
#include <random>

#include "projectUtils.h"

class HashTables
{
private:
    int numOfHashTables;
    int numberOfHyperplanes, numberOfPoints, TableSize, dim;
    std::vector<std::vector<Bucket>> hash_tables;
    std::vector<std::vector<int>> ri; // r=(-100,100)
    std::vector<std::vector<double>> t;
    std::vector<std::vector<std::vector<double>>> v;

public:
    HashTables(int L, int numberOfHyperplanes, int numberOfPoints, int dimension, int tableSize); // Constructor
    void InsertPoint(PointPtr point);
    void InsertCurve(PointPtr curve);
    int HashFunc(PointPtr point, int hashtableId);
    void PrintHashTables();
    kNeighboursPtr find_k_nearest_neighbours(PointPtr queryPoint, int k_neighbours);
    std::vector<PointPtr> range_search(PointPtr queryPoint, double range, std::vector<std::string> *foundPoints = NULL);
};

#endif
