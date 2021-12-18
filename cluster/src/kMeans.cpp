#include "kMeans.h"

// Input: Vector of all input Points
// Output: Vector of k centroid points
std::vector<PointPtr> k_means(std::vector<PointPtr> inputPoints, clusterInputData *CLData) //int numOfCentroidPoints, int dimension)
{

    std::vector<PointPtr> centroidPoints;

    int numOfInputPoints = CLData->numberOfInputPoints; //inputPoints.size();

    PointPtr currPoint;

    // choose first centroid point at random
    currPoint = inputPoints[floor(uniformDistributionGenerator(0.0, numOfInputPoints * 1.0))];
    numOfInputPoints--;
    centroidPoints.push_back(currPoint);
    std::remove(inputPoints.begin(), inputPoints.end(), currPoint);

    for (int t = 1; t < CLData->number_of_clusters; t++)
    { // first point already chosen, so t=1

        std::vector<double> D;
        D.resize(numOfInputPoints);

        for (int i = 0; i < numOfInputPoints; i++)
        {
            D[i] = min_dist_from_centroid(inputPoints[i], centroidPoints, CLData->dimension);
        }

        currPoint = inputPoints[choose_point(inputPoints, D, numOfInputPoints)];
        centroidPoints.push_back(currPoint);
        std::remove(inputPoints.begin(), inputPoints.end(), currPoint);
        numOfInputPoints--;
    }
    return centroidPoints;
}

double min_dist_from_centroid(PointPtr point, std::vector<PointPtr> centroidPoints, int dimension /*, method=l2_dist */)
{

    double min_dist = INT32_MAX;
    double cur_dist = 0.0;

    for (int i = 0; i < centroidPoints.size(); i++)
    {
        cur_dist = euclideanDistance(point, centroidPoints[i], dimension);
        if (cur_dist < min_dist)
            min_dist = cur_dist;
    }

    return min_dist;
}

void get_centroid_point(std::vector<PointPtr> inputPoints, std::vector<PointPtr> centroidPoints, PointPtr point)
{
    centroidPoints.push_back(point);
    std::remove(inputPoints.begin(), inputPoints.end(), point);
}

/* Returns index of point chosen as new centroid */
int choose_point(std::vector<PointPtr> inputPoints, std::vector<double> D, int numOfInputPoints)
{

    std::vector<double> P; // probability
    P.resize(numOfInputPoints + 1);
    P[0] = 0;

    for (int i = 1; i <= numOfInputPoints; i++)
    {
        P[i] = P[i - 1] + D[i - 1] * D[i - 1];
    }

    double x = uniformDistributionGenerator(0.0, P[numOfInputPoints]);

    // lower_bound finds in logn(binary search ?) the element after which x would go
    // distance returns dist between index of elem begin() and index of elem found
    // so we get the index of elem found
    return std::distance(P.begin(), std::lower_bound(P.begin(), P.end(), x));
}