#ifndef _PROJECT_UTILS_H_
#define _PROJECT_UTILS_H_

#include <string>
#include <vector>
#include <limits.h>

#define BIGM (4294967291)
#define W (5)
#define DELTA (3)

typedef struct PointStruct *PointPtr;
typedef struct BucketStruct *BucketPtr;
typedef struct NeighbourStruct *NeighbourPtr;
typedef struct kNeighboursStruct *kNeighboursPtr;

typedef struct PointStruct
{
    std::string id;
    std::vector<double> coords;
} Point;

typedef Point Curve;
typedef PointPtr CurvePtr;
typedef struct BucketStruct
{
    std::vector<PointPtr> points;
    std::vector<long int> ID;

} Bucket;

typedef struct NeighbourStruct
{
    PointPtr point;
    double dist;
} Neighbour;

typedef struct kNeighboursStruct
{
    std::vector<NeighbourPtr> neighbours;
    int size; // number of requested (k) nearest neighbours
} kNeighbours;

struct BY_ID
{
    // Source: https://stackoverflow.com/questions/2999135/how-can-i-sort-the-vector-elements-using-members-as-the-key-in-c
    bool operator()(PointPtr const &a, PointPtr const &b) const
    {
        return a->id < b->id;
    }
};

struct BY_ID_INT
{
    // Source: https://stackoverflow.com/questions/2999135/how-can-i-sort-the-vector-elements-using-members-as-the-key-in-c
    bool operator()(PointPtr const &a, PointPtr const &b) const
    {
        return stoi(a->id) < stoi(b->id);
    }
};
std::vector<std::string> get_lines(std::string fileName);
int get_points(std::vector<std::string> linesVector, std::vector<PointPtr> *pointVector);
int get_curves(std::vector<std::string> linesVector, std::vector<CurvePtr> *curvesVector);
// std::vector<PointPtr> *convert_points(int dimension, const std::vector<PointPtr> *pointVector);
std::vector<CurvePtr> *convert_points(int dimension, const std::vector<PointPtr> *pointVector);
void sort_neighbours(kNeighboursPtr k_nearest_neighbours, int k_neighbours);
void sort_points(std::vector<PointPtr> *Data);
void sort_points_str(std::vector<std::string> *Data);
int notAlreadyExists(kNeighboursPtr k_nearest_neighbours, std::string pointID);
kNeighboursPtr find_k_true_neighbours(PointPtr queryPoint, int k_neighbours, std::vector<PointPtr> inputPoints, int dimension);
std::string checkRerun();
CurvePtr snap_curve(const CurvePtr curve, double delta, std::vector<double> *taf, int dimension);

PointPtr snap_point(const PointPtr point, int delta, int dimension);
PointPtr concat_point(const PointPtr point, int dimension);
void remove_dup_points(PointPtr point, int dimension);
void pad_curve(CurvePtr curve, int dim);

#endif