#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <algorithm> //std::find

#include "lsh_frechet_cont.h"
#include "lshUtils.h"
#include "mathUtils.h"
#include "projectUtils.h"

FrechetContinuousHashTables::FrechetContinuousHashTables(int L, int numberOfHyperplanes, int numberOfPoints, int dimension, int tableSize) // Constructor

{
    this->numOfHashTables = L;
    this->numberOfHyperplanes = numberOfHyperplanes;
    this->numberOfPoints = numberOfPoints;
    this->TableSize = tableSize;

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

void FrechetContinuousHashTables::FrContInsertPoint(PointPtr point)
{
    for (int i = 0; i < this->numOfHashTables; i++)
    {
        PointPtr concated_point = new PointStruct;
        concated_point->coords.resize(this->dim * 2);
        std::vector<PointPtr> _curve, snapped_curve;

        pointToCurve(point, &_curve, this->dim);
        filter_curve(&_curve, this->dim, EPSILON);
        snap_curve_cont(&snapped_curve, &_curve, this->_delta, this->dim);
        minimaximize_curve_cont(&snapped_curve, this->dim);
        pad_curve_new(&snapped_curve, this->dim);
        curveToPoint(concated_point, &snapped_curve, this->dim);

        long id = FrContHashFunc(concated_point, i);
        int j = euclideanModulo(id, this->TableSize);

        this->hash_tables[i][j].ID.push_back(FrContHashFunc(concated_point, i));
        this->hash_tables[i][j].points.push_back(point);

        deleteCrv(&snapped_curve);
        deleteCrv(&_curve);
        snapped_curve.clear();
        _curve.clear();

        delete concated_point;
    }
}

int FrechetContinuousHashTables::FrContHashFunc(PointPtr point, int hashtableId)
{
    long h;
    long hri = 0;
    for (int i = 0; i < numberOfHyperplanes; i++)
    {
        // hri += this->ri[hashtableId][i] * floor((inner_product(point->coords.begin(), point->coords.end(), this->v[hashtableId][i].begin(), 0) + this->t[hashtableId][i]) / W);
        h = floor((inner_product(point->coords.begin(), point->coords.end(), this->v[hashtableId][i].begin(), 0) + this->t[hashtableId][i]) / W);

        hri += avoidOverFlowModulo(this->ri[hashtableId][i], h, BIGM, '*');
    }
    return euclideanModulo(hri, BIGM);
}

kNeighboursPtr FrechetContinuousHashTables::FrCont_find_k_nearest_neighbours(PointPtr queryPoint, int k_neighbours)
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

    for (int i = 0; i < this->numOfHashTables; i++) // for i from 1 to L do
    {
        PointPtr concated_point = concat_point(queryPoint, this->dim);
        crv _curve, snapped_curve;
        pointToCurve(concated_point, &_curve, this->dim);
        filter_curve(&_curve, this->dim, EPSILON);
        snap_curve_cont(&snapped_curve, &_curve, this->_delta, this->dim);
        minimaximize_curve_cont(&snapped_curve, this->dim);
        pad_curve_new(&snapped_curve, this->dim);
        curveToPoint(concated_point, &snapped_curve, this->dim);

        long queryID = this->FrContHashFunc(concated_point, i);
        int g = euclideanModulo(queryID, this->TableSize);
        for (int j = 0; j < this->hash_tables[i][g].points.size(); j++) // for each item p in bucket gi(q) do
        {

            if (this->hash_tables[i][g].ID[j] == queryID && notAlreadyExists(returnData, this->hash_tables[i][g].points[j]->id)) // if p,q actually belong in same bucket
            {

                currNeighbour->point = this->hash_tables[i][g].points[j];
                currNeighbour->dist = ContinuousFrechetDistance(originalQueryPoint, currNeighbour->point, this->dim);
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

    if (returnData->size < k_neighbours * numOfHashTables)
    { // rerun without ID check if haven't found enough neighbors

        for (int i = 0; i < this->numOfHashTables; i++) // for i from 1 to L do
        {
            PointPtr concated_point = concat_point(queryPoint, this->dim);
            crv _curve, snapped_curve;
            pointToCurve(concated_point, &_curve, this->dim);
            filter_curve(&_curve, this->dim, EPSILON);
            snap_curve_cont(&snapped_curve, &_curve, this->_delta, this->dim);
            minimaximize_curve_cont(&snapped_curve, this->dim);
            pad_curve_new(&snapped_curve, this->dim);
            curveToPoint(concated_point, &snapped_curve, this->dim);

            long queryID = this->FrContHashFunc(concated_point, i);
            int g = euclideanModulo(queryID, this->TableSize);

            for (int j = 0; j < this->hash_tables[i][g].points.size(); j++) // for each item p in bucket gi(q) do
            {

                if (this->hash_tables[i][g].ID[j] != queryID && notAlreadyExists(returnData, this->hash_tables[i][g].points[j]->id)) // if p,q actually belong in same bucket
                {

                    currNeighbour->point = this->hash_tables[i][g].points[j];
                    currNeighbour->dist = ContinuousFrechetDistance(originalQueryPoint, currNeighbour->point, this->dim);
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

double ContinuousFrechetDistance(PointPtr p, PointPtr q, int dimension)
{
    Curve *fp = convertToFredCurve(p, dimension);
    Curve *fq = convertToFredCurve(q, dimension);
    struct Frechet::Continuous::Distance dist = Frechet::Continuous::distance(*fp, *fq);
    delete fp;
    delete fq;
    return dist.value;
    // turn point p to Fred Curve
    // turn point q to Fred Curve
    // call Fred Cont Dist
    // return Fred Cont Dist
}

Curve *convertToFredCurve(PointPtr p, int dim)
{
    Points fp(1);
    for (int i = 0; i < dim; i++)
    {
        Point t(1);
        t.set(0, p->coords[i * 2 + 1]);
        fp.add(t);
    }

    Curve *curve = new Curve(fp);
    return curve;
}