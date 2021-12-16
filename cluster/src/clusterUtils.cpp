#include <fstream>
#include <algorithm>

#include "clusterUtils.h"
#include "lsh_frechet_dsc.h"
#include "methods.h"

double silhouette_calculator(PointPtr point, std::vector<Cluster> *clusters, int dimension)
{

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

std::vector<int> get_2_closest_clusters(PointPtr point, std::vector<Cluster> *clusters, int dimension)
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
    std::vector<PointPtr>(*centroidPoints);
    for (int i = 0; i < (*clusters).size(); i++)
    {
        (*centroidPoints).push_back((*clusters)[i].centroidPoint);
    }

    for (int i = 0; i < (*centroidPoints).size(); i++)
    {

        cur_dist = euclideanDistance(point, (*centroidPoints)[i], dimension);
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

double calculateChanges(std::vector<PointPtr> *centroids, std::vector<Cluster> *clusters, std::vector<PointPtr> *newCentroids, int dimension)
{
    int numOfClusters = clusters->size();

    double change = 0.0;

    if (newCentroids == NULL)
        newCentroids = new std::vector<PointPtr>;
    else
        newCentroids->clear();

    newCentroids->resize(numOfClusters);

    for (int i = 0; i < numOfClusters; i++)
    {
        (*newCentroids)[i] = new Point(*((*centroids)[i]));
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
    std::vector<treeNodePtr> trees;
    trees.resize(clusters->size());
    std::vector<int> inputed;
    double changes = 0;
    for (int i = 0; i < clusters->size(); i++)
    {
        trees[i] = buildTree(floor(log2((*clusters)[i].points.size())));
        fillTree(trees[i], (*clusters)[i].points, &inputed);
        inputed.clear();
        (*newCentroids)[i] = findMean(trees[i]);
        changes += DFDistance((*centroids)[i], (*newCentroids)[i], dimension * 2);
    }

    return change;
}

PointPtr computeMeanCurve(PointPtr curve1, PointPtr curve2)
{
    std::vector<std::vector<double>> _c;
    double dis = DFDistance(curve1, curve2, curve1.size(), &_c);
    std::vector<std::vector<int>> traversal;
    computeOptimalTraversal(&_c, &traversal, curve1->coords.size() / 2);
    PointPtr retPoint = New PointStruct;
    for (int i = 0; i < curve1->coords.size() / 2; i++)
    {
        double ii = (curve1->coords[traversal[i][0]] + curve2->coords[traversal[i][0]]) / 2;
        double jj = (curve1->coords[traversal[i][1]] + curve2->coords[traversal[i][1]]) / 2;
        retPoint->coords.push_back(ii);
        retPoint->coords.push_back(jj);
    }
    return retPoint;
}

void computeOptimalTraversal(std::vector<std::vector<double>> *_c, std::vector<std::vector<int>> *traversal, int dimension)
{
    std::vector<std::vector<int>> traversal;
    std::vector<int> element;
    element.resize(2);
    int index_p = dimension - 1, index_q = dimension - 1;
    element[0] = index_p, element[1] = index_q;
    traversal->push_back(element);
    while (index_q != 0 && index_p != 0)
    {
        int minIdx = minIdx(_c[index_p - 1, index_q], _c[index_p, index_q - 1], _c[index_p - 1, index_q - 1]); // Fix
        if (minIdx == 0)
        {
            element[0] = --index_p, element[1] = index_q;
            traversal->push_back(element);
        }
        else if (minIdx == 0)
        {
            element[0] = index_p, element[1] = --index_q;
            traversal->push_back(element);
        }
        else
        {
            element[0] = --index_p, element[1] = --index_q;
            traversal->push_back(element);
        }
    }
    std::reverse(traversal->begin(), traversal->end());
}

PointPtr findMean(treeNodePtr treeNode)
{
    if (treeNode->curve != NULL)
        return treeNode->curve;
    PointPtr meanLeft = findMean(treeNode->leftChld);
    PointPtr meanRight = findMean(treeNode->rightChld);
    PointPtr mean = computeMeanCurve(treeNode->leftChld, treeNode->rightChld);
    return mean;
}

void fillTree(treeNodePtr treeNode, std::vector<PointPtr> c_points, std::vector<int> *inp)
{
    if (treeNode->leftChld == NULL && treeNode->rightChld == NULL)
    {
        int index;
        do
            index = (int)ceil(rand() % c_points.size());
        while (std::find(inp->begin(), inp->end(), index) != inp->end());

        treeNode->curve = c_points[index];
        inp->push_back(index);
    }
    else
    {
        treeNode->curve = NULL;
        fillTree(treeNode->leftChld, c_points, inp);
        fillTree(treeNode->rightChld, c_points, inp);
    }
}

int getInputData(int argc, char **argv, inputData *CLData)
{
    std::vector<std::string> found;
    CLData = new inputData;

    // Initializing with default values, which may change depending on the config file's content
    CLData->number_of_vector_hash_tables = DEF_VECTOR_HASH_TABLES;
    CLData->number_of_vector_hash_functions = DEF_VECTOR_HASH_FUNCTIONS;
    CLData->max_number_M_hypercube = DEF_MAX_NUM_M_CUBE;
    CLData->number_of_hypercube_dimensions = DEF_NUM_CUBE_DIM;
    CLData->number_of_probes = DEF_PROBES;
    CLData->complete = false;

    for (int i = 0; i < argc; i++)
    {

        if (std::string(argv[i]) == "-i")
        {
            CLData->inputFileName = std::string(argv[i + 1]);
            std::cout << CLData->inputFileName << std::endl;
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
        else if (std::string(argv[i]) == "-m")
        {
            CLData->methodName = argv[i + 1];
            std::cout << CLData->methodName << std::endl;
            found.push_back("m");
        }
    }

    found.push_back(" ");

    if (std::find(found.begin(), found.end(), "inputFile") == found.end()) // if not found inputFile
    {
        std::cerr << "Input file name not given! Please try again using -i <input file>" << std::endl;
        return EXIT_FAIL_INPUT_ERR;
    }
    if (std::find(found.begin(), found.end(), "configFile") == found.end()) // if not found inputFile
    {
        std::cerr << "Config file name not given! Please try again using -c <config file>" << std::endl;
        return EXIT_FAIL_CONFIG_ERR;
    }
    if (std::find(found.begin(), found.end(), "outputFile") == found.end()) // if not found inputFile
    {
        std::cerr << "Output file name not given! Please try again using -o <output file>" << std::endl;
        return EXIT_FAIL_OUTPUT_ERR;
    }
    if (std::find(found.begin(), found.end(), "m") == found.end()) // if not found inputFile
    {
        std::cerr << "Method name not given! Please try again using -m <Classic OR LSH OR Hypercube>" << std::endl;
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
        else
        {
            std::cerr << "Invalid method name! Please try again using -m <Classic OR LSH OR Hypercube>" << std::endl;
            return EXIT_FAIL_METHOD_ERR;
        }
    }
    found.clear();
    std::ifstream configFile(CLData->configFileName);
    if (!configFile.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << CLData->configFileName << "'"
                  << std::endl;
        return EXIT_FAIL_CONFIG_ERR;
    }
    std::cout << "Reading config file " << CLData->configFileName << "..." << std::endl;

    std::string line;

    while (getline(configFile, line))
    {
        // inputLines.push_back(line);
        std::string word = "";
        std::string parameter = "";
        int value = -1;
        for (char x : line)
        {
            if (x == ' ')
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
    configFile.close();

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

int execCluster(inputData *CLData, std::vector<Cluster> *clusters, std::vector<PointPtr> *inputPoints, std::vector<PointPtr> *centroidPoints)
{
    std::vector<PointPtr> tempCentroidPoints;
    bool flag = false;
    if (CLData->method == CLASSIC_METHOD)
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
            change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension);

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
    else if (CLData->method == LSH_METHOD)
    {
        HashTables HashTablesObject(CLData->number_of_vector_hash_tables, CLData->number_of_vector_hash_functions, CLData->numberOfInputPoints, CLData->dimension, CLData->numberOfInputPoints / 8);

        for (int i = 0; i < CLData->numberOfInputPoints; i++)
            HashTablesObject.HashTables::InsertPoint(((*inputPoints))[i]);
        double change = INT32_MAX * 1.0;
        int count = 0;
        while (change >= TOL && count < 30)
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
            change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension);

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
    else if (CLData->method == HYPERCUBE_METHOD)
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
            change = calculateChanges(centroidPoints, clusters, &tempCentroidPoints, CLData->dimension);

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
    return EXIT_SUCCESS;
}

double evalSilhouette(inputData *CLData, std::vector<Cluster> *clusters)
{
    double totalSilhouette = 0.0;
    for (int i = 0; i < CLData->number_of_clusters; i++)
    { // for each cluster
        double silhouetteSum = 0.0;
        for (int j = 0; j < (*clusters)[i].size; j++)
        { // for each point in cluster
            silhouetteSum += silhouette_calculator((*clusters)[i].points[j], clusters, CLData->dimension);
        }
        (*clusters)[i].silhouette = silhouetteSum / (double)((*clusters)[i].size); // saves average
        totalSilhouette += silhouetteSum;
    }
    totalSilhouette /= CLData->numberOfInputPoints;

    return totalSilhouette;
}

int writeToOutput(inputData *CLData, std::vector<Cluster> *clusters, std::vector<PointPtr> *centroidPoints, double totalSilhouette, int tCluster)
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
    outputFile << "Silhouette: [";

    for (int i = 0; i < CLData->number_of_clusters; i++)
        outputFile << (*clusters)[i].silhouette << ", ";
    outputFile << totalSilhouette << "]" << std::endl;

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

void deleteData(std::vector<PointPtr> *inputPoints, inputData *CLData, std::vector<PointPtr> *centroidPoints)
{
    for (int i = 0; i < CLData->numberOfInputPoints; i++)
        delete (*inputPoints)[i];

    for (int i = 0; CLData->number_of_clusters; i++)
        delete (*centroidPoints)[i];

    delete CLData;
}