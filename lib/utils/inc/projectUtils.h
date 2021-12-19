#ifndef _PROJECT_UTILS_H_
#define _PROJECT_UTILS_H_

#include <string>
#include <vector>
#include <limits.h>

#define BIGM (4294967291)
#define W (30)
#define DELTA (1.73)
#define EPSILON (1)

typedef struct Point__Struct *PointPtr;
typedef struct BucketStruct *BucketPtr;
typedef struct NeighbourStruct *NeighbourPtr;
typedef struct kNeighboursStruct *kNeighboursPtr;
typedef std::vector<PointPtr> crv;
typedef std::vector<PointPtr> *crvPtr;

typedef struct Point__Struct
{
    std::string id;
    std::vector<double> coords;
} PointStruct;

// typedef Point Curve;
// typedef PointPtr CurvePtr;
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

typedef struct treeNode treeNode;
typedef struct treeNode *treeNodePtr;
typedef struct treeNode
{
    PointPtr curve;
    treeNodePtr rightChld, leftChld;
} treeNode;

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
int get_curves(std::vector<std::string> linesVector, std::vector<PointPtr> *curvesVector);
// std::vector<PointPtr> *convert_points(int dimension, const std::vector<PointPtr> *pointVector);
std::vector<PointPtr> *convert_points(int dimension, const std::vector<PointPtr> *pointVector);
void sort_neighbours(kNeighboursPtr k_nearest_neighbours, int k_neighbours);
void sort_points(std::vector<PointPtr> *Data);
void sort_points_str(std::vector<std::string> *Data);
int notAlreadyExists(kNeighboursPtr k_nearest_neighbours, std::string pointID);
kNeighboursPtr find_k_true_neighbours(PointPtr queryPoint, int k_neighbours, std::vector<PointPtr> inputPoints, int dimension, int useEuclDist = 1);
kNeighboursPtr find_k_true_neighbours_dfd(PointPtr queryPoint, int k_neighbours, std::vector<PointPtr> inputPoints, int dim);
std::string checkRerun();
crvPtr snap_curve(crvPtr snapped_curve, const crvPtr curve, double delta, std::vector<double> *taf, int dimension);

/// Aii ///
PointPtr snap_point(const PointPtr point, int delta, int dimension);

PointPtr concat_point(const PointPtr point, int dimension);
void remove_dup_points(crvPtr curve, int dimension);
void pad_curve(PointPtr curve, int dim);
void pad_curve_new(crvPtr curve, int dim);

/// Aiii ///
void filter_curve(crvPtr curve, int dimension, double epsilon);
crvPtr snap_curve_cont(crvPtr snapped_curve, crvPtr curve, double delta, int dimension);
void minimaximize_curve_cont(crvPtr _curve, int dimension);

void pointToCurve(const PointPtr _p, crvPtr _c, int dimension);
void curveToPoint(PointPtr _p, crvPtr _c, int dimension);

void deleteCrv(crvPtr _curve);

treeNodePtr buildTree(double height);

#endif