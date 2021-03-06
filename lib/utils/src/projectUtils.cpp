#include <fstream>
#include <unistd.h>

#include "../inc/projectUtils.h"
#include "../inc/mathUtils.h"

std::vector<std::string> get_lines(std::string fileName)
{
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << fileName << "'"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::vector<std::string> inputLines;
    std::string line;
    while (std::getline(file, line))
    {
        inputLines.push_back(line);
    }
    file.close();
    return inputLines;
}

int get_points(std::vector<std::string> linesVector, std::vector<PointPtr> *pointsVector)
{
    int dimension;
    for (int i = 0; i < linesVector.size(); i++)
    {
        // separate std::string by Tabs
        // pick every element from 2nd to std::endl
        // read point coordinates
        PointPtr currPoint = new PointStruct;
        std::string word = "";
        dimension = 0;
        for (char x : linesVector[i])
        {
            if (x == '\t')
            {
                if (dimension)
                    currPoint->coords.push_back(atof(word.c_str()));
                else
                    currPoint->id = word;
                word = "";

                dimension++;
            }
            else
            {
                word = word + x;
            }
        }
        currPoint->coords.push_back(atof(word.c_str()));
        pointsVector->push_back(currPoint);
    }
    return dimension;
}

PointPtr snap_point(const PointPtr point, int delta, int dimension)
{
    int t;
    PointPtr snapped_point = new PointStruct;
    snapped_point->id = point->id;
    for (double x : point->coords)
    {
        t = uniformDistributionGenerator(0.0, delta);
        snapped_point->coords.push_back(floor(abs(x - t) / delta + 0.5) * delta + t);
    }

    return snapped_point;
}

void sort_neighbours(kNeighboursPtr k_nearest_neighbours, int k_neighbours) // sort distance in a vector of k distances
{
    int k = k_neighbours; // number of neighbours
    NeighbourPtr tempNeighbour;

    for (int i = k - 1; i > 0; i--)
    {
        if (k_nearest_neighbours->neighbours[i]->dist < k_nearest_neighbours->neighbours[i - 1]->dist)
        {
            tempNeighbour = k_nearest_neighbours->neighbours[i - 1];
            k_nearest_neighbours->neighbours[i - 1] = k_nearest_neighbours->neighbours[i];
            k_nearest_neighbours->neighbours[i] = tempNeighbour;
        }
    }
}

void sort_points(std::vector<PointPtr> *Data) // sort distance in a vector of k distances
{
    int k = Data->size();
    PointPtr tempPoint;
    for (int i = k - 1; i > 0; i--)
    {
        if ((*Data)[i]->id < (*Data)[i - 1]->id)
        {
            tempPoint = (*Data)[i];
            (*Data)[i] = (*Data)[i - 1];
            (*Data)[i - 1] = tempPoint;
        }
    }
}

void sort_points_str(std::vector<std::string> *Data)
{
    int k = Data->size();
    std::string tempPoint;
    for (int i = k - 1; i > 0; i--)
    {
        if ((*Data)[i] < (*Data)[i - 1])
        {
            tempPoint = (*Data)[i];
            (*Data)[i] = (*Data)[i - 1];
            (*Data)[i - 1] = tempPoint;
        }
    }
}

int notAlreadyExists(kNeighboursPtr k_nearest_neighbours, std::string pointID)
{

    for (int i = 0; i < k_nearest_neighbours->size; i++)
        if (k_nearest_neighbours->neighbours[i]->point->id == pointID)
            return 0;
    return 1;
}

kNeighboursPtr find_k_true_neighbours(PointPtr queryPoint, int k_neighbours, std::vector<PointPtr> inputPoints, int dim, int useEuclDist)
{
    if (!useEuclDist)
        return find_k_true_neighbours_dfd(concat_point(queryPoint, dim), k_neighbours, inputPoints, dim);

    NeighbourPtr currNeighbour = new Neighbour;
    kNeighboursPtr returnData = new kNeighbours;

    returnData->neighbours.resize(k_neighbours);
    returnData->size = k_neighbours;

    for (int i = 0; i < k_neighbours; i++)
    {
        returnData->neighbours[i] = new Neighbour;
        returnData->neighbours[i]->point = NULL;
        returnData->neighbours[i]->dist = INT32_MAX; // initialize distance with a very big value
    }

    for (int i = 0; i < inputPoints.size(); i++)
    {
        currNeighbour->point = inputPoints[i];
        currNeighbour->dist = euclideanDistance(queryPoint, currNeighbour->point, dim);

        if (currNeighbour->dist < returnData->neighbours[k_neighbours - 1]->dist && currNeighbour->dist > 0)
        {
            // if (returnData->size < k_neighbours)
            //     returnData->size++;
            returnData->neighbours[k_neighbours - 1]->point = currNeighbour->point;
            returnData->neighbours[k_neighbours - 1]->dist = currNeighbour->dist;

            if (k_neighbours > 1)
                sort_neighbours(returnData, k_neighbours);
        }
    }

    delete currNeighbour;
    return returnData;
}

kNeighboursPtr find_k_true_neighbours_dfd(PointPtr queryPoint, int k_neighbours, std::vector<PointPtr> inputPoints, int dim)
{
    int numOfInputPoints = inputPoints.size();
    std::vector<PointPtr> *inputPoints_2d = new std::vector<PointPtr>;
    inputPoints_2d->resize(numOfInputPoints);
    for (int i = 0; i < numOfInputPoints; i++)
    {
        (*inputPoints_2d)[i] = concat_point(inputPoints[i], dim);
    }

    NeighbourPtr currNeighbour = new Neighbour;
    kNeighboursPtr returnData = new kNeighbours;

    returnData->neighbours.resize(k_neighbours);
    returnData->size = k_neighbours;

    for (int i = 0; i < k_neighbours; i++)
    {
        returnData->neighbours[i] = new Neighbour;
        returnData->neighbours[i]->point = NULL;
        returnData->neighbours[i]->dist = INT32_MAX; // initialize distance with a very big value
    }

    for (int i = 0; i < numOfInputPoints; i++)
    {
        currNeighbour->point = (*inputPoints_2d)[i];
        currNeighbour->dist = DFDistance(queryPoint, currNeighbour->point, dim);

        if (currNeighbour->dist < returnData->neighbours[k_neighbours - 1]->dist && currNeighbour->dist > 0)
        {
            // if (returnData->size < k_neighbours)
            //     returnData->size++;
            returnData->neighbours[k_neighbours - 1]->point = currNeighbour->point;
            returnData->neighbours[k_neighbours - 1]->dist = currNeighbour->dist;

            if (k_neighbours > 1)
                sort_neighbours(returnData, k_neighbours);
        }
    }

    delete currNeighbour;
    return returnData;
}

std::string checkRerun()
{
    std::cout << "Rerun Program?..." << std::endl
              << "======Options======" << std::endl
              << "CONT to rerun" << std::endl
              << "TERM to terminate" << std::endl
              << "===================" << std::endl;
    std::string option;
    char ch;
    while (true)
    {
        option = "";
        while ((ch = getchar()) != '\n')
            option += ch;
        if (option == "CONT" || option == "TERM" || option == "cont" || option == "term")
            break;
    }
    return option;
}

// Save input file lines as a discretized curve,
// projected from 2D space to a single vector of size 2*d.
int get_curves(std::vector<std::string> linesVector, std::vector<PointPtr> *pointsVector)
{
    int dimension;
    for (int i = 0; i < linesVector.size(); i++)
    {
        // separate std::string by Tabs
        // pick every element from 2nd to std::endl
        // read point coordinates
        PointPtr currPoint = new PointStruct;
        std::string word = "";
        dimension = 0;
        for (char x : linesVector[i])
        {
            if (x == '\t')
            {
                if (dimension)
                    currPoint->coords.push_back(atof(word.c_str()));
                else
                    currPoint->id = word;
                word = "";

                dimension++;
            }
            else
            {
                word = word + x;
            }
        }
        currPoint->coords.push_back(atof(word.c_str()));
        pointsVector->push_back(currPoint);
    }

    return dimension;
}

// // Convert euclidean vector point to discretized curve,
// // projected from 2D space to a single vector of size 2*d.
// std::vector<PointPtr> *convert_points(int dimension, const std::vector<PointPtr> *pointsVector)
// {
//     std::vector<PointPtr> *converted_points = new std::vector<PointPtr>;
//     for (int i = 0; i < pointsVector->size(); i++)
//     {
//         PointPtr point = new Point;
//         point->id = (*pointsVector)[i]->id;
//         for (int j = 0; j < dimension; j++)
//         {
//             point->coords.push_back((*pointsVector)[i]->coords[j]);
//             point->coords.push_back(j);
//         }
//         converted_points->push_back(point);
//     }
//     return converted_points;
// }

// Convert euclidean vector point to discretized curve,
// projected from 2D space to a single vector of size 2*d.
std::vector<PointPtr> *convert_points(int dimension, const std::vector<PointPtr> *point_vector)
{
    std::vector<PointPtr> *curves = new std::vector<PointPtr>;
    for (int i = 0; i < point_vector->size(); i++)
    {
        PointPtr curve = new PointStruct;
        curve->id = (*point_vector)[i]->id;
        for (int j = 0; j < dimension; j++)
        {
            curve->coords.push_back((*point_vector)[i]->coords[j]); // xi
            curve->coords.push_back(j);                             // yi
        }
        curves->push_back(curve);
    }
    return curves;
}

crvPtr snap_curve(crvPtr snapped_curve, const crvPtr curve, double delta, std::vector<double> *taf, int dimension)
{
    snapped_curve->resize(dimension);
    for (int i = 0; i < dimension; i++)
    {
        (*snapped_curve)[i] = new PointStruct;
        (*snapped_curve)[i]->id = (*curve)[i]->id;
        (*snapped_curve)[i]->coords.resize(2);
        (*snapped_curve)[i]->coords[0] = (floor(abs((*curve)[i]->coords[0] - (*taf)[0]) / delta + 0.5) * delta + (*taf)[0]);
        (*snapped_curve)[i]->coords[1] = (floor(abs((*curve)[i]->coords[1] - (*taf)[1]) / delta + 0.5) * delta + (*taf)[1]);
    }

    return snapped_curve;
}

crvPtr snap_curve_cont(crvPtr snapped_curve, crvPtr curve, double delta, int dimension)
{
    snapped_curve->resize(dimension);
    for (int i = 0; i < dimension; i++)
    {
        (*snapped_curve)[i] = new PointStruct;
        (*snapped_curve)[i]->id = (*curve)[i]->id;
        (*snapped_curve)[i]->coords.resize(2);
        (*snapped_curve)[i]->coords[0] = floor((*curve)[i]->coords[0] / delta) * delta;
        (*snapped_curve)[i]->coords[1] = floor((*curve)[i]->coords[1] / delta) * delta;
    }

    return snapped_curve;
}

PointPtr concat_point(const PointPtr point, int dimension)
{
    PointPtr concatd_point = new PointStruct;
    concatd_point->id = point->id;
    concatd_point->coords.resize(2 * dimension);
    for (int i = 0; i < dimension; i++)
    {
        concatd_point->coords[i * 2] = i;
        concatd_point->coords[(i * 2) + 1] = point->coords[i];
    }

    return concatd_point;
}

void remove_dup_points(crvPtr curve, int dimension)
{
    std::vector<int> removedIndex;
    double prev[2];
    prev[0] = (*curve)[0]->coords[0];
    prev[1] = (*curve)[0]->coords[1];
    for (int i = 1; i < dimension; i++)
    {
        if ((*curve)[i]->coords[0] == prev[0] && (*curve)[i]->coords[1] == prev[1])
        {
            removedIndex.push_back(i);
        }
        else
        {
            prev[0] = (*curve)[i]->coords[0];
            prev[1] = (*curve)[i]->coords[1];
        }
    }
    int i = 0;
    // for (int index : removedIndex)
    // {
    //     // delete (*curve)[index];
    //     curve->erase(curve->begin() + index);
    // }

    for (int i = 0; i < removedIndex.size(); i++)
    {
        curve->erase(curve->begin() + removedIndex[i] - i);
    }
}

void pad_curve(PointPtr curve, int dim)
{
    int curveSize = curve->coords.size();
    for (int i = curveSize; i < dim; i++)
    {

        curve->coords.push_back(INT_MAX - 10);
    }
}

void pad_curve_new(crvPtr curve, int dim)
{
    int curveSize = curve->size();
    for (int i = curveSize - 1; i < dim; i++)
    {
        PointPtr _point = new PointStruct;
        _point->coords.resize(2);
        _point->coords[0] = 10000.0;
        _point->coords[1] = 10000.0;
        curve->push_back(_point);
    }
}

//////////////// Aiii ////////////////

void filter_curve(crvPtr curve, int dimension, double epsilon)
{
    // Calculate points that need to be removed
    for (int i = 1, index = 1; i < dimension - 1; i++)
    {
        if (euclideanDistance((*curve)[index], (*curve)[index - 1], 2) <= epsilon && euclideanDistance((*curve)[index], (*curve)[index + 1], 2) <= epsilon)
        {
            delete (*curve)[index];
            curve->erase(curve->begin() + index);
        }
        else
            index++;
    }
}

void minimaximize_curve_cont(crvPtr _curve, int dimension)
{
    std::vector<int> removedIndex;
    for (int i = 1; i < dimension - 1; i++)
    {
        if ((*_curve)[i]->coords[1] > (*_curve)[i - 1]->coords[1] && (*_curve)[i]->coords[1] < (*_curve)[i + 1]->coords[1])
        {
            removedIndex.push_back(i);
        }
    }

    for (int i = 0; i < removedIndex.size(); i++)
    {
        _curve->erase(_curve->begin() + removedIndex[i] - i);
    }
}

void pointToCurve(const PointPtr _p, crvPtr _c, int dimension)
{
    if (_p == NULL)
    {
        std::cerr << "Structs have not been initialized!\n Cannot make conversion!" << std::endl;
        return;
    }

    _c->resize(dimension);
    for (int i = 0; i < dimension; i++)
    {
        (*_c)[i] = new PointStruct;
        (*_c)[i]->id = _p->id;
        (*_c)[i]->coords.resize(2);
        (*_c)[i]->coords[0] = _p->coords[i * 2];
        (*_c)[i]->coords[1] = _p->coords[(i * 2) + 1];
    }
}

void curveToPoint(PointPtr _p, crvPtr _c, int dimension)
{
    for (int i = 0; i < dimension; i++)
    {
        _p->coords[i * 2] = (*_c)[i]->coords[0];
        _p->coords[(i * 2) + 1] = (*_c)[i]->coords[1];
    }
}

void deleteCrv(crvPtr _curve)
{
    for (int i = 0; i < _curve->size(); i++)
    {
        delete (*_curve)[i];
    }
}

treeNodePtr buildTree(int height)
{
    if (height <= 0)
    {
        return NULL;
    }
    treeNodePtr retTree = new treeNode;
    retTree->curve = NULL;
    retTree->rightChld = buildTree(height - 1);
    retTree->leftChld = buildTree(height - 1);
    return retTree;
}