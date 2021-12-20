#include <fstream>
#include <algorithm>
#include <ctime>

#include "clusterUtils.h"
#include "lsh_frechet_dsc.h"
#include "methods.h"

int counter = 0;

void frechet_method(FrechetDiscreteHashTables *HashTablesObject, std::vector<PointPtr> *centroids, std::vector<Cluster> *clusters, const std::vector<PointPtr> *inputPoints, clusterInputData *CLData, int numOfInputPoints)
{

    std::vector<PointPtr> foundPoints;
    std::vector<std::string> foundPointIDs;
    std::vector<std::vector<std::string>> foundPointIDsPerCluster;

    foundPointIDsPerCluster.resize(CLData->number_of_clusters);

    int inputPointsSize = CLData->numberOfInputPoints;
    double currRadius = minFrDistBetweenCentroids(centroids, CLData->number_of_clusters, CLData->dimension) / 2;
    std::vector<std::vector<PointPtr>> clusterPoints;
    std::vector<PointPtr> duplicates;
    clusterPoints.resize(CLData->number_of_clusters);
    int initialInputPoints = inputPointsSize;
    double initialRadius = currRadius;
    int prevNumOfFound = 0;
    int currNumOfFound = 0;
    int numOfFound = 0;
    // std::cout << "initialInput points" << initialInputPoints << std::endl;
    // std::cout << "numOfFound" << numOfFound << std::endl;
    // std::cout << "currad" << currRadius << std::endl;
    // std::cout << "InitialRadius" << initialRadius << std::endl;
    // std::cout << "currNumOfFound" << currNumOfFound << std::endl;
    // std::cout << "prevNumOfFound" << prevNumOfFound << std::endl;
    int loopCounter = 0;
    while (initialInputPoints - numOfFound >= initialInputPoints / 2 && currRadius < initialRadius * 100 && (currNumOfFound >= prevNumOfFound || currNumOfFound > 1) && loopCounter < 5)
    {
        loopCounter++;
        prevNumOfFound = currNumOfFound;
        currNumOfFound = 0;
        for (int c = 0; c < CLData->number_of_clusters; c++)
        {
            clusterPoints[c] = HashTablesObject->range_search((*centroids)[c], currRadius, &(foundPointIDsPerCluster[c]));
        }
        std::vector<std::string> tempArray;
        for (int i = 0; i < CLData->number_of_clusters; i++)
        {
            std::merge(foundPointIDsPerCluster[i].begin(), foundPointIDsPerCluster[i].end(), tempArray.begin(), tempArray.end(), std::back_inserter(foundPointIDs));
            tempArray.clear();
            for (auto currPoint : foundPointIDs)
            {
                tempArray.push_back(currPoint);
            }
        }
        for (int i = 0; i < CLData->number_of_clusters; i++)
        {
            foundPointIDsPerCluster[i].clear();
            for (auto currPoint : foundPointIDs)
            {
                foundPointIDsPerCluster[i].push_back(currPoint);
            }
        }

        duplicates = find_duplicates(clusterPoints, CLData->number_of_clusters);
        for (auto currPoint : duplicates)
        {
            std::vector<int> CentroidsToBeErased;
            std::vector<int> position;
            for (int c = 0; c < CLData->number_of_clusters; c++)
            {
                for (int p = 0; p < clusterPoints[c].size(); p++)
                {

                    if (currPoint->id == clusterPoints[c][p]->id)
                    {
                        CentroidsToBeErased.push_back(c);
                        position.push_back(p);
                    }
                }
            }
            double minDist = INT_MAX;
            int minCentroid;
            double currDist = 0.0;
            int currCentroid;

            for (int i = 0; i < CentroidsToBeErased.size(); i++)
            {
                currCentroid = i;
                currDist = DFDistance((*centroids)[currCentroid], currPoint, CLData->dimension);
                if (currDist < minDist)
                {
                    minDist = currDist;
                    minCentroid = currCentroid;
                }
            }
            for (int i = 0; i < CentroidsToBeErased.size(); i++)
            {
                if (i != minCentroid)
                {
                    clusterPoints[CentroidsToBeErased[i]].erase(clusterPoints[CentroidsToBeErased[i]].begin() + position[i]);
                }
            }
        }
        for (int c = 0; c < CLData->number_of_clusters; c++)
        {
            for (int pointIndex = 0; pointIndex < clusterPoints[c].size(); pointIndex++)
            {
                (*clusters)[c].points.push_back(clusterPoints[c][pointIndex]);
                (*clusters)[c].size++;
                currNumOfFound++;
            }
            clusterPoints[c].clear();
        }
        // std::sort(foundPointIDs.begin(), foundPointIDs.end());
        currRadius *= 2;
        numOfFound += currNumOfFound;
    }
    int index = 0;

    for (auto currPoint : (*inputPoints))
    {
        bool found = false;
        for (int c = 0; c < CLData->number_of_clusters; c++)
        {
            // search for inputPoint in every cluster
            if (find((*clusters)[c].points.begin(), (*clusters)[c].points.end(), currPoint) != (*clusters)[c].points.end())
            {
                found = true;
                break;
            }
        }
        if (!found)
        {
            // if not in any clusters, find closest manually
            index = lloyd_method(centroids, currPoint, CLData->dimension);
            (*clusters)[index].points.push_back(currPoint);
            (*clusters)[index].size++;
        }
    }
    foundPoints.clear();
    foundPointIDs.clear();
    for (int c = 0; c < CLData->number_of_clusters; c++)
        foundPointIDsPerCluster[c].clear();
    foundPointIDsPerCluster.clear();
}

double minDistBetweenCentroids(std::vector<PointPtr> *centroidPoints, int numOfCentroids, int dimension)
{
    double minDist = INT_MAX;
    double currDist = 0.0;

    for (int i = 0; i < numOfCentroids; i++)
    {
        for (int j = 0; j < i; j++)
        {
            currDist = euclideanDistance((*centroidPoints)[i], (*centroidPoints)[j], dimension);
            if (currDist < minDist)
                minDist = currDist;
        }
    }
    return minDist;
}

double minFrDistBetweenCentroids(std::vector<PointPtr> *centroidPoints, int numOfCentroids, int dimension)
{
    double minDist = INT_MAX;
    double currDist = 0.0;

    for (int i = 0; i < numOfCentroids; i++)
    {
        for (int j = 0; j < i; j++)
        {
            currDist = DFDistance((*centroidPoints)[i], (*centroidPoints)[j], dimension);
            if (currDist < minDist)
                minDist = currDist;
        }
    }
    return minDist;
}

std::vector<PointPtr> find_duplicates(std::vector<std::vector<PointPtr>> clusterPoints, int numOfClusters)
{
    std::vector<PointPtr> mainArray;
    std::vector<PointPtr> tempArray;
    std::vector<PointPtr> duplPoints;
    for (int i = 0; i < numOfClusters; i++)
    {
        mainArray.clear();
        std::merge((clusterPoints)[i].begin(), (clusterPoints)[i].end(), tempArray.begin(), tempArray.end(), std::back_inserter(mainArray), BY_ID());
        tempArray.clear();
        for (auto currPoint : mainArray)
        {
            tempArray.push_back(currPoint);
        }
    }

    if (!mainArray.empty())
    {
        PointPtr currPoint = mainArray[0];
        int i = 1;
        while (i < mainArray.size())
        {
            if (currPoint->id == mainArray[i]->id)
            {
                duplPoints.push_back(currPoint);
                i++;
                while (i < mainArray.size())
                {
                    if (currPoint->id != mainArray[i]->id)
                        break;
                    i++;
                }
                continue;
            }
            currPoint = mainArray[i];
            i++;
        }
    }
    return duplPoints;
}

double silhouette_calculator(PointPtr point, std::vector<Cluster> *clusters, int dimension, int method)
{

    if (method == UPDATE_FRECHET)
        return silhouette_calculator_frechet(point, clusters, dimension);

    // find 2 closest (*clusters)
    std::vector<int> closestClusters = get_2_closest_clusters(point, clusters, dimension);

    // find ai
    ClusterPtr cluster = &((*clusters)[closestClusters[0]]);
    double distanceSum = 0;
    for (int i = 0; i < cluster->size; i++)
    {
        distanceSum += euclideanDistance(point, cluster->points[i], dimension);
    }
    double a = distanceSum / (cluster->size * 1.0); // -1 beacause point belongs in cluster

    // find bi
    cluster = &((*clusters)[closestClusters[1]]);
    distanceSum = 0;
    for (int i = 0; i < cluster->size; i++)
    {
        distanceSum += euclideanDistance(point, cluster->points[i], dimension);
    }
    double b = distanceSum / (cluster->size * 1.0); // -1 beacause point belongs in cluster

    return (b - a) / (std::max(a, b));
}

double silhouette_calculator_frechet(PointPtr point, std::vector<Cluster> *clusters, int dimension)
{
    // find 2 closest (*clusters)
    std::vector<int> closestClusters = get_2_closest_clusters(point, clusters, dimension, 0);
    // find ai
    ClusterPtr cluster = &((*clusters)[closestClusters[0]]);
    double distanceSum = 0;
    for (int i = 0; i < cluster->size; i++)
    {
        distanceSum += DFDistance(point, cluster->points[i], dimension);
    }
    double a = distanceSum / (cluster->size * 1.0); // -1 beacause point belongs in cluster

    // find bi
    cluster = &((*clusters)[closestClusters[1]]);
    distanceSum = 0;
    for (int i = 0; i < cluster->size; i++)
    {
        distanceSum += DFDistance(point, cluster->points[i], dimension);
    }
    double b = distanceSum / (cluster->size * 1.0); // -1 beacause point belongs in cluster

    return (b - a) / (std::max(a, b));
}

std::vector<int> get_2_closest_clusters(PointPtr point, std::vector<Cluster> *clusters, int dimension, int useEuclDist)
{

    std::vector<int> closestClusters;
    closestClusters.resize(2);
    for (int i = 0; i < 2; i++)
        closestClusters[i] = 0;

    std::vector<double> closestDist;
    closestDist.resize(2);
    for (int i = 0; i < 2; i++)
        closestDist[i] = INT32_MAX;

    double min_dist = INT32_MAX;
    double cur_dist = 0.0;
    // get all centroid points ( we could add this as a parameter)
    std::vector<PointPtr> centroidPoints;
    for (int i = 0; i < (*clusters).size(); i++)
    {
        centroidPoints.push_back((*clusters)[i].centroidPoint);
    }

    for (int i = 0; i < centroidPoints.size(); i++)
    {
        if (useEuclDist)
            cur_dist = euclideanDistance(point, centroidPoints[i], dimension);
        else
            cur_dist = DFDistance(point, centroidPoints[i], dimension);
        if (cur_dist < closestDist[1])
        {
            closestDist[1] = cur_dist;
            closestClusters[1] = i;
        }
        // sort 1&2
        if (closestDist[1] < closestDist[0])
        {
            closestDist[1] = closestDist[0];
            closestDist[0] = cur_dist; // cur_dist holds temp
            closestClusters[1] = closestClusters[0];
            closestClusters[0] = i; // i holds cur_cluster
        }
    }

    return closestClusters;
}

double calculateChanges(std::vector<PointPtr> *centroids, std::vector<Cluster> *clusters, std::vector<PointPtr> *newCentroids, int dimension, int method)
{
    if (method == UPDATE_FRECHET)
    {
        return calculateChangesCurve(centroids, clusters, newCentroids, dimension);
    }
    int numOfClusters = clusters->size();

    double change = 0.0;

    if (newCentroids == NULL)
        newCentroids = new std::vector<PointPtr>;
    else
        newCentroids->clear();

    newCentroids->resize(numOfClusters);

    for (int i = 0; i < numOfClusters; i++)
    {
        (*newCentroids)[i] = new PointStruct(*((*centroids)[i]));
        int numOfPoints = (*clusters)[i].points.size();
        for (int j = 0; j < dimension; j++)
        {
            (*newCentroids)[i]->coords[j] = 0.0;
            for (int k = 0; k < numOfPoints; k++)
            {
                (*newCentroids)[i]->coords[j] += (*clusters)[i].points[k]->coords[j];
            }
            (*newCentroids)[i]->coords[j] /= numOfPoints;
        }

        // cluster->centroidPoint = centroid;
        change += euclideanDistance((*centroids)[i], (*newCentroids)[i], dimension);
    }

    return change;
}

double calculateChangesCurve(std::vector<PointPtr> *centroids, std::vector<Cluster> *clusters, std::vector<PointPtr> *newCentroids, int dimension)
{
    if (newCentroids == NULL)
        newCentroids = new std::vector<PointPtr>;
    else
        newCentroids->clear();

    newCentroids->resize(clusters->size());

    std::vector<treeNodePtr> trees;
    trees.resize(clusters->size());
    std::vector<int> inputed;
    double changes = 0;
    for (int i = 0; i < clusters->size(); i++)
    {
        if ((*clusters)[i].points.size() != 0)
        {
            std::vector<PointPtr> safeClusterPoints;
            safeClusterPoints.resize((*clusters)[i].points.size());
            for (int j = 0; j < (*clusters)[i].points.size(); j++)
            {
                safeClusterPoints[j] = new PointStruct;
                safeClusterPoints[j]->coords.resize(dimension * 2);
                safeClusterPoints[j]->id = (*clusters)[i].points[j]->id;
                for (int k = 0; k < dimension * 2; k++)
                    safeClusterPoints[j]->coords[k] = (*clusters)[i].points[j]->coords[k];
            }
            int height = ceil(log2(safeClusterPoints.size()));
            trees[i] = buildTree(height);
            fillTree(trees[i], safeClusterPoints, &inputed);
            counter = 0;
            inputed.clear();
            (*newCentroids)[i] = findMean(trees[i]);
            changes += DFDistance((*centroids)[i], (*newCentroids)[i], dimension);
        }
        else
        {
            (*newCentroids)[i] = new PointStruct;
            (*newCentroids)[i]->id = (*centroids)[i]->id;
            for (int jj = 0; jj < dimension * 2; jj++)
            {
                (*newCentroids)[i]->coords.push_back((*centroids)[i]->coords[jj]);
            }
        }
    }

    return changes;
}

PointPtr computeMeanCurve(PointPtr curve1, PointPtr curve2)
{
    if (curve2 == NULL)
    {
        return curve1;
    }
    std::vector<std::vector<double>> *_c = new std::vector<std::vector<double>>;

    double dis = DFDistance(curve1, curve2, curve1->coords.size() / 2, _c);
    std::vector<std::vector<int>> traversal;
    computeOptimalTraversal(_c, &traversal, curve1->coords.size() / 2);
    PointPtr retPoint = new PointStruct;
    retPoint->id = "0";
    for (int i = 0; i < curve1->coords.size() / 2; i++)
    {
        double ii = (curve1->coords[traversal[i][0]] + curve2->coords[traversal[i][0]]) / 2;
        double jj = (curve1->coords[traversal[i][1]] + curve2->coords[traversal[i][1]]) / 2;
        retPoint->coords.push_back(ii);
        retPoint->coords.push_back(jj);
    }
    if (_c != NULL)
        delete _c;
    return retPoint;
}

void computeOptimalTraversal(std::vector<std::vector<double>> *_c, std::vector<std::vector<int>> *traversal, int dimension)
{
    // std::vector<std::vector<int>> traversal;
    std::vector<int> element;
    element.resize(2);
    int index_p = dimension - 1, index_q = dimension - 1;
    element[0] = index_p, element[1] = index_q;
    traversal->push_back(element);
    while (index_q != 0 && index_p != 0)
    {
        int minIdx = minIndx((*_c)[index_p - 1][index_q], (*_c)[index_p][index_q - 1], (*_c)[index_p - 1][index_q - 1]); // Fix

        if (minIdx == 0)
        {
            element.clear();
            element.resize(2);
            element[0] = --index_p, element[1] = index_q;
            traversal->push_back(element);
        }
        else if (minIdx == 0)
        {
            element.clear();
            element.resize(2);
            element[0] = index_p, element[1] = --index_q;
            traversal->push_back(element);
        }
        else
        {
            element.clear();
            element.resize(2);
            element[0] = --index_p, element[1] = --index_q;
            traversal->push_back(element);
        }
    }
    std::reverse(traversal->begin(), traversal->end());
}

PointPtr findMean(treeNodePtr treeNode)
{
    PointPtr meanLeft, meanRight;
    if (treeNode->leftChld == NULL && treeNode->rightChld == NULL)
    {
        return treeNode->curve;
    }
    else
    {
        PointPtr meanLeft = findMean(treeNode->leftChld);
        if (treeNode->rightChld != NULL)
            PointPtr meanRight = findMean(treeNode->rightChld);
        else
            PointPtr meanRight = NULL;
        PointPtr mean = computeMeanCurve(meanLeft, meanRight);
        return mean;
    }
}

void fillTree(treeNodePtr treeNode, std::vector<PointPtr> c_points, std::vector<int> *inp)
{
    if (treeNode == NULL)
    {
        return;
    }
    if (treeNode->leftChld == NULL && treeNode->rightChld == NULL)
    {
        treeNode->curve = c_points[counter++];
    }
    else
    {
        treeNode->curve = NULL;
        fillTree(treeNode->leftChld, c_points, inp);
        fillTree(treeNode->rightChld, c_points, inp);
    }
}

int getCLInputData(int argc, char **argv, clusterInputData *CLData)
{
    std::vector<std::string> found;

    // Initializing with default values, which may change depending on the config file's content
    CLData->number_of_clusters = 3;
    CLData->number_of_vector_hash_tables = DEF_VECTOR_HASH_TABLES;
    CLData->number_of_vector_hash_functions = DEF_VECTOR_HASH_FUNCTIONS;
    CLData->max_number_M_hypercube = DEF_MAX_NUM_M_CUBE;
    CLData->number_of_hypercube_dimensions = DEF_NUM_CUBE_DIM;
    CLData->number_of_probes = DEF_PROBES;
    CLData->complete = false;
    CLData->silhouette = false;

    for (int i = 0; i < argc; i++)
    {

        if (std::string(argv[i]) == "-i")
        {
            CLData->inputFileName = std::string(argv[i + 1]);
            std::cout << "inputFileName: " << CLData->inputFileName << std::endl;
            found.push_back("inputFile");
        }
        else if (std::string(argv[i]) == "-c")
        {
            CLData->configFileName = std::string(argv[i + 1]);
            std::cout << CLData->configFileName << std::endl;
            found.push_back("configFile");
        }
        else if (std::string(argv[i]) == "-o")
        {
            CLData->outputFileName = std::string(argv[i + 1]);
            std::cout << CLData->outputFileName << std::endl;
            found.push_back("outputFile");
        }
        else if (std::string(argv[i]) == "-complete")
        {
            CLData->complete = true;
            std::cout << "-complete" << std::endl;
            found.push_back("complete");
        }
        else if (std::string(argv[i]) == "-silhouette")
        {
            CLData->silhouette = true;
            std::cout << "-silhouette" << std::endl;
            found.push_back("silhouette");
        }
        else if (std::string(argv[i]) == "-assignment")
        {
            CLData->methodName = argv[i + 1];
            std::cout << CLData->methodName << std::endl;
            found.push_back("assignment");
        }
        else if (std::string(argv[i]) == "-update")
        {
            CLData->updateName = argv[i + 1];
            std::cout << CLData->updateName << std::endl;
            found.push_back("update");
        }
    }

    found.push_back(" ");

    if (std::find(found.begin(), found.end(), "inputFile") == found.end()) // if not found inputFile
    {
        std::cerr << "Input file name not given! Please try again using -i <input file>" << std::endl;
        return EXIT_FAIL_INPUT_ERR;
    }
    if (std::find(found.begin(), found.end(), "configFile") == found.end()) // if not found configFile
    {
        std::cerr << "Config file name not given! Please try again using -c <config file>" << std::endl;
        return EXIT_FAIL_CONFIG_ERR;
    }
    if (std::find(found.begin(), found.end(), "outputFile") == found.end()) // if not found outputFile
    {
        std::cerr << "Output file name not given! Please try again using -o <output file>" << std::endl;
        return EXIT_FAIL_OUTPUT_ERR;
    }
    if (std::find(found.begin(), found.end(), "assignment") == found.end()) // if not found assignment
    {
        std::cerr << "Assignment not given! Please try again using -assignment <Classic OR LSH OR Hypercube OR LSH_Frechet>" << std::endl;
        return EXIT_FAIL_METHOD_ERR;
    }
    if (std::find(found.begin(), found.end(), "update") == found.end()) // if not found update
    {
        std::cerr << "Update method name not given! Please try again using -update <Mean Frechet or Mean Vector>" << std::endl;
        return EXIT_FAIL_METHOD_ERR;
    }
    else
    {
        if (CLData->methodName == "Classic")
        {
            CLData->method = CLASSIC_METHOD;
        }
        else if (CLData->methodName == "LSH")
        {
            CLData->method = LSH_METHOD;
        }
        else if (CLData->methodName == "Hypercube")
        {
            CLData->method = HYPERCUBE_METHOD;
        }
        else if (CLData->methodName == "LSH_Frechet")
        {
            CLData->method = FRECHET_D_METHOD;
        }
        else
        {
            std::cerr << "Invalid method name! Please try again using -m <Classic OR LSH OR Hypercube OR LSH_Frechet>" << std::endl;
            return EXIT_FAIL_METHOD_ERR;
        }

        if (CLData->updateName == "Mean_Frechet")
        {
            CLData->update = UPDATE_FRECHET;
        }
        else if (CLData->updateName == "Mean_Vector")
        {
            CLData->update = UPDATE_VECTOR;
        }
    }
    found.clear();

    std::vector<std::string> conflines = get_lines(CLData->configFileName);

    for (std::string line : conflines)
    {
        line.push_back('\n');
        std::string word = "";
        std::string parameter = "";
        int value = -1;
        for (char x : line)
        {
            if (x == ' ' || x == '\n')
            {
                if (parameter.empty())
                {
                    word.pop_back();
                    parameter = word;
                    found.push_back(word);
                }
                else if (value == -1)
                {
                    if (!is_number(word))
                    {
                        std::cerr << "Parameter [" << parameter << "]: Value '" << word << "' is not an integer" << std::endl;
                        return EXIT_FAIL_CONFIG_ERR;
                    }

                    value = stoi(word);
                    break;
                }
                word = "";
            }
            else
            {
                word = word + x;
            }
        }

        if (parameter == "number_of_clusters")
        {
            CLData->number_of_clusters = value;
        }
        else if (parameter == "number_of_vector_hash_tables")
        {
            CLData->number_of_vector_hash_tables = value;
        }
        else if (parameter == "number_of_vector_hash_functions")
        {
            CLData->number_of_vector_hash_functions = value;
        }
        else if (parameter == "max_number_M_hypercube")
        {
            CLData->max_number_M_hypercube = value;
        }
        else if (parameter == "number_of_hypercube_dimensions")
        {
            CLData->number_of_hypercube_dimensions = value;
        }
        else if (parameter == "number_of_probes")
        {
            CLData->number_of_probes = value;
        }
    }

    std::cout << "Config file parsed" << std::endl;

    found.push_back(" ");

    if (std::find(found.begin(), found.end(), "number_of_clusters") == found.end()) // if not found number_of_clusters
    {
        std::cerr << "Config file must contain parameter [number_of_clusters]."
                  << std::endl
                  << "Please include 'number_of_clusters: <int>' in the config file..."
                  << std::endl;
        return EXIT_FAIL_CONFIG_ERR;
    }
    return EXIT_SUCCESS;
}

int execCluster(clusterInputData *CLData, std::vector<Cluster> *clusters, std::vector<PointPtr> *inputPoints, std::vector<PointPtr> *inputPoints_2d, std::vector<PointPtr> *centroidPoints)
{
    std::vector<PointPtr> tempCentroidPoints;
    bool flag = false;
    if (CLData->method == CLASSIC_METHOD)
    {
        if (CLData->update == UPDATE_VECTOR)
        {
            double change = INT32_MAX * 1.0;
            int count = 0;
            while (change >= TOL / 1000 && count < 50)
            {
                if (flag)
                {
                    for (int i = 0; i < CLData->number_of_clusters; i++)
                    {
                        delete (*centroidPoints)[i];
                        (*centroidPoints)[i] = tempCentroidPoints[i];
                        (*clusters)[i].centroidPoint = tempCentroidPoints[i];
                    }
                    tempCentroidPoints.clear();
                }
                else
                    flag = true;
                for (int c = 0; c < CLData->number_of_clusters; c++)
                {
                    (*clusters)[c].points.clear();
                    (*clusters)[c].size = 0;
                }
                int index = 0;
                for (int i = 0; i < CLData->numberOfInputPoints; i++)
                {
                    index = lloyd_method(centroidPoints, (*inputPoints)[i], CLData->dimension);
                    (*clusters)[index].points.push_back((*inputPoints)[i]);
                    (*clusters)[index].size++;
                }
                change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension, CLData->update);

                std::cout << "Change " << change << "," << count << std::endl;
                int totalPoints = 0;
                for (int i = 0; i < CLData->number_of_clusters; i++)
                {
                    totalPoints += (*clusters)[i].size;
                }
                std::cout << "Total points: " << totalPoints << std::endl;
                count++;
            }
        }
        else if (CLData->update == UPDATE_FRECHET)
        {
            double change = INT32_MAX * 1.0;
            int count = 0;
            while (change >= 1 && count < 10)
            {
                if (flag)
                {
                    for (int i = 0; i < CLData->number_of_clusters; i++)
                    {
                        delete (*centroidPoints)[i];
                        (*centroidPoints)[i] = tempCentroidPoints[i];
                        (*clusters)[i].centroidPoint = tempCentroidPoints[i];
                    }
                    tempCentroidPoints.clear();
                }
                else
                    flag = true;
                for (int c = 0; c < CLData->number_of_clusters; c++)
                {
                    (*clusters)[c].points.clear();
                    (*clusters)[c].size = 0;
                }
                int index = 0;
                for (int i = 0; i < CLData->numberOfInputPoints; i++)
                {
                    index = lloyd_method_DFD(centroidPoints, (*inputPoints_2d)[i], CLData->dimension);
                    (*clusters)[index].points.push_back((*inputPoints_2d)[i]);
                    (*clusters)[index].size++;
                }
                change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension, CLData->update);

                std::cout << "Change " << change << "," << count << std::endl;
                int totalPoints = 0;
                for (int i = 0; i < CLData->number_of_clusters; i++)
                {
                    totalPoints += (*clusters)[i].size;
                }
                std::cout << "Total points: " << totalPoints << std::endl;
                count++;
            }
        }
    }
    else if (CLData->method == LSH_METHOD)
    {
        if (CLData->update == UPDATE_VECTOR)
        {
            HashTables HashTablesObject(CLData->number_of_vector_hash_tables, CLData->number_of_vector_hash_functions, CLData->numberOfInputPoints, CLData->dimension, CLData->numberOfInputPoints / 8);

            for (int i = 0; i < CLData->numberOfInputPoints; i++)
                HashTablesObject.HashTables::InsertPoint(((*inputPoints))[i]);
            double change = INT32_MAX * 1.0;
            int count = 0;
            while (change >= TOL && count < 10)
            {
                if (flag)
                {
                    for (int i = 0; i < CLData->number_of_clusters; i++)
                    {
                        delete (*centroidPoints)[i];
                        (*centroidPoints)[i] = tempCentroidPoints[i];
                        (*clusters)[i].centroidPoint = tempCentroidPoints[i];
                    }
                    tempCentroidPoints.clear();
                }
                else
                    flag = true;
                for (int c = 0; c < CLData->number_of_clusters; c++)
                {
                    (*clusters)[c].points.clear();
                    (*clusters)[c].size = 0;
                }
                lsh_method(&HashTablesObject, centroidPoints, clusters, inputPoints, CLData, CLData->numberOfInputPoints);
                change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension, CLData->update);

                std::cout << "Change " << change << "," << count << std::endl;
                int totalPoints = 0;
                for (int i = 0; i < CLData->number_of_clusters; i++)
                {
                    totalPoints += (*clusters)[i].size;
                }
                std::cout << "Total points: " << totalPoints << std::endl;
                count++;
            }
        }
        else
        {
            std::cerr << "Cannot use Classic LSH assignment with time curve update" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if (CLData->method == HYPERCUBE_METHOD)
    {
        if (CLData->update == UPDATE_VECTOR)
        {
            HChashTable HypercubeObject(CLData->dimension, CLData->number_of_hypercube_dimensions, CLData->number_of_probes, CLData->max_number_M_hypercube);

            for (int i = 0; i < CLData->numberOfInputPoints; i++)
                HypercubeObject.HChashTable::InsertPoint((*inputPoints)[i]);
            double change = INT32_MAX * 1.0;
            int count = 0;
            while (change > TOL && count < 30)
            {
                if (flag)
                {
                    for (int i = 0; i < CLData->number_of_clusters; i++)
                    {
                        delete (*centroidPoints)[i];
                        (*centroidPoints)[i] = tempCentroidPoints[i];
                        (*clusters)[i].centroidPoint = tempCentroidPoints[i];
                    }
                    tempCentroidPoints.clear();
                }
                else
                    flag = true;
                for (int c = 0; c < CLData->number_of_clusters; c++)
                {
                    (*clusters)[c].points.clear();
                    (*clusters)[c].size = 0;
                }
                hyperCube_method(&HypercubeObject, centroidPoints, clusters, inputPoints, CLData, CLData->numberOfInputPoints);
                change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension, CLData->update);

                std::cout << "Change " << change << "," << count << std::endl;
                int totalPoints = 0;
                for (int i = 0; i < CLData->number_of_clusters; i++)
                {
                    totalPoints += (*clusters)[i].size;
                }
                std::cout << "Total points: " << totalPoints << std::endl;
                count++;
            }
        }
        else
        {
            std::cerr << "Cannot use Classic HyperCube assignment with time curve update" << std::endl;
            return EXIT_FAILURE;
        }
    }
    else if (CLData->method == FRECHET_D_METHOD)
    {
        if (CLData->update == UPDATE_FRECHET)
        {
            FrechetDiscreteHashTables HashTablesObject(CLData->number_of_vector_hash_tables, CLData->number_of_vector_hash_functions, CLData->numberOfInputPoints, CLData->dimension, CLData->numberOfInputPoints / 8);

            for (int i = 0; i < CLData->numberOfInputPoints; i++)
                HashTablesObject.FrechetDiscreteHashTables::FrDscInsertPoint(((*inputPoints_2d))[i]);
            double change = INT32_MAX * 1.0;
            int count = 0;
            while (change >= 1 && count < 10)
            {
                if (flag)
                {
                    for (int i = 0; i < CLData->number_of_clusters; i++)
                    {
                        delete (*centroidPoints)[i];
                        (*centroidPoints)[i] = tempCentroidPoints[i];
                        (*clusters)[i].centroidPoint = tempCentroidPoints[i];
                    }
                    tempCentroidPoints.clear();
                }
                else
                    flag = true;
                for (int c = 0; c < CLData->number_of_clusters; c++)
                {
                    (*clusters)[c].points.clear();
                    (*clusters)[c].size = 0;
                }
                frechet_method(&HashTablesObject, centroidPoints, clusters, inputPoints_2d, CLData, CLData->numberOfInputPoints);
                change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension, CLData->update);

                std::cout << "Change " << change << "," << count << std::endl;
                int totalPoints = 0;
                for (int i = 0; i < CLData->number_of_clusters; i++)
                {
                    totalPoints += (*clusters)[i].size;
                }
                std::cout << "Total points: " << totalPoints << std::endl;
                count++;
            }
        }
        else
        {
            std::cerr << "Cannot use LSH Frechet assignment with vector update" << std::endl;
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}

double evalSilhouette(clusterInputData *CLData, std::vector<Cluster> *clusters)
{
    double totalSilhouette = 0.0;
    for (int i = 0; i < CLData->number_of_clusters; i++)
    { // for each cluster
        double silhouetteSum = 0.0;
        for (int j = 0; j < (*clusters)[i].size; j++)
        { // for each point in cluster
            silhouetteSum += silhouette_calculator((*clusters)[i].points[j], clusters, CLData->dimension, CLData->update);
        }
        (*clusters)[i].silhouette = silhouetteSum / (double)((*clusters)[i].size); // saves average
        totalSilhouette += silhouetteSum;
    }
    totalSilhouette /= CLData->numberOfInputPoints;
    return totalSilhouette;
}

int writeToOutput(clusterInputData *CLData, std::vector<Cluster> *clusters, std::vector<PointPtr> *centroidPoints, double totalSilhouette, int tCluster)
{
    std::ofstream outputFile(CLData->outputFileName);
    if (!outputFile.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << CLData->outputFileName << "'"
                  << std::endl;
        return EXIT_FAIL_OUTPUT_ERR;
    }

    outputFile << "Algorithm: ";
    if (CLData->method == CLASSIC_METHOD)
        outputFile << "Lloyds";
    else if (CLData->method == LSH_METHOD)
        outputFile << "Range Search LSH";
    else if (CLData->method == HYPERCUBE_METHOD)
        outputFile << "Range Search Hypercube";
    outputFile << std::endl;

    for (int i = 0; i < CLData->number_of_clusters; i++)
    {
        outputFile << "CLUSTER-"
                   << i + 1 << "{size: " << (*clusters)[i].size
                   << ", centroid: ";

        for (int j = 0; j < CLData->dimension; j++)
            outputFile << (*centroidPoints)[i]->coords[j] << " ";
        outputFile << "}" << std::endl
                   << std::endl;
    }
    outputFile << "clustering_time: " << (double)(tCluster / 1000)
               << "s" << std::endl;

    if (CLData->silhouette)
    {
        outputFile << "Silhouette: [";

        for (int i = 0; i < CLData->number_of_clusters; i++)
            outputFile << (*clusters)[i].silhouette << ", ";
        outputFile << totalSilhouette << "]" << std::endl;
    }

    if (CLData->complete)
    {
        outputFile << std::endl;
        for (int i = 0; i < CLData->number_of_clusters; i++)
        {
            std::sort((*clusters)[i].points.begin(), (*clusters)[i].points.end(), BY_ID_INT());
            outputFile << std::endl
                       << "CLUSTER-"
                       << i + 1 << " {"
                       << (*clusters)[i].centroidPoint->id;

            for (int j = 0; j < (*clusters)[i].size; j++)
                outputFile << ", " << (*clusters)[i].points[j]->id;

            outputFile << "}" << std::endl;
        }
    }

    outputFile.close();

    return EXIT_SUCCESS;
}

void deleteData(std::vector<PointPtr> *centroidPoints, std::vector<PointPtr> *inputPoints, std::vector<PointPtr> *inputPoints_2d, clusterInputData *CLData)
{
    for (int i = 0; i < CLData->numberOfInputPoints; i++)
    {
        if (inputPoints_2d != NULL)
            if ((*inputPoints_2d)[i] != NULL)
                delete (*inputPoints_2d)[i];
    }

    for (int i = 0; i < CLData->number_of_clusters; i++)
        if ((*centroidPoints)[i] != NULL)
            delete (*centroidPoints)[i];

    delete CLData;
}