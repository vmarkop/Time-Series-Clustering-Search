#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <algorithm> //std::find

#include "lsh_frechet_dsc.h"

FrechetDiscreteHashTables::FrechetDiscreteHashTables(int L, int numberOfHyperplanes, int numberOfPoints, int dimension, int tableSize) // Constructor

{
    this->numOfHashTables = L;
    this->numberOfHyperplanes = numberOfHyperplanes;
    this->numberOfPoints = numberOfPoints;
    this->TableSize = tableSize;
    std::cout << "TABLA" << L << std::endl;

    this->dim = dimension;
    this->hash_tables.resize(L);
    this->t.resize(L);
    this->ri.resize(L);
    this->v.resize(L);
    this->_taf.resize(L);
    this->_delta = DELTA;

    for (int i = 0; i < L; i++)
    {
        this->_taf[i].resize(2);
        this->_taf[i][0] = uniformDistributionGenerator(0.0, this->_delta);
        this->_taf[i][1] = uniformDistributionGenerator(0.0, this->_delta);

        this->hash_tables[i].resize(this->TableSize);

        std::cout << "TABLA" << L << std::endl;
        this->hash_tables[i][0].ID.push_back(1);
        this->hash_tables[i][0].ID.clear();

        std::cout << "TABLA" << L << std::endl;
        this->t[i].resize(numberOfHyperplanes);
        this->ri[i].resize(numberOfHyperplanes);
        this->v[i].resize(numberOfHyperplanes);

        for (int j = 0; j < numberOfHyperplanes; j++)
        {

            this->t[i][j] = uniformDistributionGenerator(0.0, W * 1.0);
            this->ri[i][j] = rand() % 2000 - 1000;
            this->v[i][j].resize(this->dim * 2);
            for (int l = 0; l < this->dim * 2; l++)
                this->v[i][j][l] = normalDistributionGenerator(0.0, 1.0);
        }
    }
}

void FrechetDiscreteHashTables::FrDscInsertPoint(PointPtr point)
{
    // (1,x),(2,y),(3,z)
    PointPtr original_point = concat_point(point, this->dim);
    PointPtr concated_point = concat_point(point, this->dim);
    for (int i = 0; i < this->numOfHashTables; i++)
    {
        std::vector<PointPtr> _curve, snapped_curve;

        pointToCurve(concated_point, &_curve, this->dim);
        snap_curve(&snapped_curve, &_curve, this->_delta, &(this->_taf[i]), this->dim);
        remove_dup_points(&snapped_curve, this->dim);
        pad_curve_new(&snapped_curve, this->dim);
        curveToPoint(concated_point, &snapped_curve, this->dim);

        int id = FrDscHashFunc(concated_point, i);
        int j = euclideanModulo(id, this->TableSize);
        deleteCrv(&snapped_curve);
        deleteCrv(&_curve);
        snapped_curve.clear();
        _curve.clear();
        this->hash_tables[i][j].ID.push_back(id);
        this->hash_tables[i][j].points.push_back(original_point);

        //(r1h1 + r2h2 + r3h3 + r4h4 + r5h5) % m = ((r1h1 % m) + (r2h2 % m) + (r3h3 % m) + (r4h4 % m) + (r5h5 % m)) % m
        // = = = ( ((r1%m * h1%m)) % m + ... + ((r5%m * h5%m)) % m ) % m
        //
    }
    delete concated_point;
}

int FrechetDiscreteHashTables::FrDscHashFunc(PointPtr point, int hashtableId)
{
    int h;
    int hri = 0;
    for (int i = 0; i < numberOfHyperplanes; i++)
    {
        // hri += this->ri[hashtableId][i] * floor((inner_product(point->coords.begin(), point->coords.end(), this->v[hashtableId][i].begin(), 0) + this->t[hashtableId][i]) / W);
        h = floor((inner_product(point->coords.begin(), point->coords.end(), this->v[hashtableId][i].begin(), 0) + this->t[hashtableId][i]) / W);

        hri += avoidOverFlowModulo(this->ri[hashtableId][i], h, BIGM, '*');
    }
    return euclideanModulo(hri, BIGM);
}

kNeighboursPtr FrechetDiscreteHashTables::FrDsc_find_k_nearest_neighbours(PointPtr queryPoint, int k_neighbours)
{

    PointPtr originalQueryPoint = concat_point(queryPoint, this->dim);
    // PointPtr curPoint;
    // int curDist;
    int count = 0;
    NeighbourPtr currNeighbour = new Neighbour;

    kNeighboursPtr returnData = new kNeighbours;
    returnData->neighbours.resize(k_neighbours);
    returnData->size = 0;

    for (int i = 0; i < k_neighbours; i++)
    {
        returnData->neighbours[i] = new Neighbour;
        returnData->neighbours[i]->point = NULL;
        returnData->neighbours[i]->dist = INT32_MAX; // initialize distance with a very big value
    }
    std::cout << "numOfHash = " << this->numOfHashTables << std::endl;
    for (int i = 0; i < this->numOfHashTables; i++) // for i from 1 to L do
    {
        PointPtr concated_point = concat_point(queryPoint, this->dim);
        crv _curve, snapped_curve;
        pointToCurve(concated_point, &_curve, this->dim);
        snap_curve(&snapped_curve, &_curve, this->_delta, &(this->_taf[i]), this->dim);
        remove_dup_points(&snapped_curve, this->dim);
        pad_curve_new(&snapped_curve, this->dim);
        curveToPoint(concated_point, &snapped_curve, this->dim);
        int queryID = FrDscHashFunc(concated_point, i);
        int g = euclideanModulo(queryID, this->TableSize);
        for (int j = 0; j < this->hash_tables[i][g].points.size(); j++) // for each item p in bucket gi(q) do
        {
            if (this->hash_tables[i][g].ID[j] == queryID && notAlreadyExists(returnData, this->hash_tables[i][g].points[j]->id)) // if p,q actually belong in same bucket
            {
                currNeighbour->point = this->hash_tables[i][g].points[j];
                currNeighbour->dist = DFDistance(originalQueryPoint, currNeighbour->point, this->dim);
                // if dist(q,p) < db then b <- p; db <- dist(q,p)
                if (currNeighbour->dist < returnData->neighbours[k_neighbours - 1]->dist)
                {
                    if (returnData->size < k_neighbours)
                        returnData->size++;
                    returnData->neighbours[k_neighbours - 1]->point = currNeighbour->point;
                    returnData->neighbours[k_neighbours - 1]->dist = currNeighbour->dist;
                    count++;

                    sort_neighbours(returnData, k_neighbours);
                }
            }
            if (count > 20 * numOfHashTables)
            {
                delete currNeighbour;
                return returnData;
            }
        }
    }

    if (returnData->size < k_neighbours)
    { // rerun without ID check if haven't found enough neighbors

        for (int i = 0; i < this->numOfHashTables; i++) // for i from 1 to L do
        {
            PointPtr concated_point = concat_point(queryPoint, this->dim);
            crv _curve, snapped_curve;
            pointToCurve(concated_point, &_curve, this->dim);
            snap_curve(&snapped_curve, &_curve, this->_delta, &(this->_taf[i]), this->dim);
            remove_dup_points(&snapped_curve, this->dim);
            pad_curve_new(&snapped_curve, this->dim);
            curveToPoint(concated_point, &snapped_curve, this->dim);
            int queryID = FrDscHashFunc(concated_point, i);
            int g = euclideanModulo(queryID, this->TableSize);

            for (int j = 0; j < this->hash_tables[i][g].points.size(); j++) // for each item p in bucket gi(q) do
            {

                if (this->hash_tables[i][g].ID[j] != queryID && notAlreadyExists(returnData, this->hash_tables[i][g].points[j]->id)) // if p,q actually belong in same bucket
                {

                    currNeighbour->point = this->hash_tables[i][g].points[j];
                    currNeighbour->dist = euclideanDistance(originalQueryPoint, currNeighbour->point, this->dim);
                    // if dist(q,p) < db then b <- p; db <- dist(q,p)
                    if (currNeighbour->dist < returnData->neighbours[k_neighbours - 1]->dist)
                    {
                        if (returnData->size < k_neighbours)
                            returnData->size++;
                        returnData->neighbours[k_neighbours - 1]->point = currNeighbour->point;
                        returnData->neighbours[k_neighbours - 1]->dist = currNeighbour->dist;

                        count++;

                        sort_neighbours(returnData, k_neighbours);
                    }
                }
                if (count > 10 * numOfHashTables)
                {
                    delete currNeighbour;
                    return returnData;
                }
            }
        }
    }
    delete currNeighbour;
    std::cout << "returnDataSize: " << returnData->size << std::endl;
    return returnData;
}

// std::vector<PointPtr> HashTables::range_search(PointPtr queryPoint, double range, std::vector<std::string> *foundPoints)
// {

//     bool noFoundPoints = false; // flag to know if data needs to be freed or not
//     if (foundPoints == NULL)
//     {
//         noFoundPoints = true;
//         foundPoints = new std::vector<std::string>;
//     }
//     else
//         std::sort(foundPoints->begin(), foundPoints->end());

//     NeighbourPtr currNeighbour = new Neighbour;

//     std::vector<PointPtr> returnData;

//     for (int i = 0; i < this->numOfHashTables; i++) // for i from 1 to L do
//     {
//         int queryID = this->HashFunc(queryPoint, i);
//         int g = euclideanModulo(queryID, this->TableSize);
//         for (int j = 0; j < this->hash_tables[i][g].points.size(); j++) // for each item p in bucket gi(q) do
//         {
//             bool found = binary_search(foundPoints->begin(), foundPoints->end(), this->hash_tables[i][g].points[j]->id);

//             if (!found)
//             {

//                 currNeighbour->point = this->hash_tables[i][g].points[j];
//                 currNeighbour->dist = euclideanDistance(queryPoint, currNeighbour->point, this->dim);

//                 if (currNeighbour->dist < range)
//                 {
//                     returnData.push_back(this->hash_tables[i][g].points[j]);
//                     foundPoints->push_back(this->hash_tables[i][g].points[j]->id);
//                     sort_points(&returnData);
//                     sort_points_str(foundPoints);
//                 }
//             }
//         }
//     }
//     delete currNeighbour;
//     if (noFoundPoints)
//         delete foundPoints;
//     return returnData;
// }

// void HashTables::PrintHashTables()
// {

//     for (int i = 0; i < this->numOfHashTables; i++)
//     {
//         std::cout << "Hash Table " << i << ":" << std::endl;
//         for (int j = 0; j < this->TableSize; j++)
//         {
//             std::cout << "\tBucket " << j << ":" << std::endl;
//             for (PointPtr &point : this->hash_tables[i][j].points)
//                 std::cout << "\t\tPoint with id: " << point->id << std::endl;
//         }
//     }
// }
