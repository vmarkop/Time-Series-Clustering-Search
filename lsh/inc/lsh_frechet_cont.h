#ifndef _LSH_FRECHET_CONT_
#define _LSH_FRECHET_CONT_

#include <stdlib.h>
#include <vector>
#include <string>
#include <random>

#include "projectUtils.h"

#include "../../Fred/include/curve.hpp"
#include "../../Fred/include/point.hpp"
#include "../../Fred/include/frechet.hpp"

class FrechetContinuousHashTables
{
private:
    int numOfHashTables;
    int numberOfHyperplanes, numberOfPoints, TableSize, dim;
    std::vector<std::vector<Bucket>> hash_tables;
    std::vector<std::vector<int>> ri; // r=(-100,100)
    std::vector<std::vector<double>> t;
    std::vector<std::vector<std::vector<double>>> v;
    std::vector<std::vector<double>> _taf;
    double _delta;

public:
    FrechetContinuousHashTables(int L, int numberOfHyperplanes, int numberOfPoints, int dimension, int tableSize); // Constructor
    void FrContInsertPoint(PointPtr point);
    int FrContHashFunc(PointPtr point, int hashtableId);
    // void PrintHashTables();
    kNeighboursPtr FrCont_find_k_nearest_neighbours(PointPtr queryPoint, int k_neighbours);
    // std::vector<PointPtr> range_search(PointPtr queryPoint, double range, std::vector<std::string> *foundPoints = NULL);
};

double ContinuousFrechetDistance(PointPtr p, PointPtr q, int dimension);
Curve *convertToFredCurve(PointPtr p, int dim);

#endif