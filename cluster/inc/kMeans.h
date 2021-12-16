#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

#include "clusterUtils.h"
#include "../../lib/projectUtils.h"

std::vector<PointPtr> k_means(std::vector<PointPtr> inputPoints, int numOfCentroidPoints, int dimension);
double min_dist_from_centroid(PointPtr point, std::vector<PointPtr> centroidPoints, int dimension /*, method=l2_dist */);
int choose_point(std::vector<PointPtr> inputPoints, std::vector<double> D);