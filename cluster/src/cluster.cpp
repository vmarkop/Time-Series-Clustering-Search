#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <climits>
#include <string>

#include "clusterUtils.h"
#include "kMeans.h"
#include "methods.h"

#include "mathUtils.h"
#include "projectUtils.h"

int main(int argc, char **argv)
{
    clusterInputData *CLData = new clusterInputData;
    int error;
    if ((error = getCLInputData(argc, argv, CLData)) > 0)
        return error;

    std::cout << "Getting lines" << std::endl;
    std::vector<std::string> inputLines = get_lines(CLData->inputFileName);
    std::vector<PointPtr> inputPoints;
    CLData->numberOfInputPoints = inputLines.size();
    std::cout << "Getting points" << std::endl;
    CLData->dimension = get_points(inputLines, &inputPoints);

    std::cout << "Dimension:" << CLData->dimension << std::endl;

    // Calculate vector of centroid points (1 for each cluster)
    std::cout << "Calculating centroid points..." << CLData->number_of_clusters << std::endl;
    std::vector<PointPtr> centroidPoints;
    centroidPoints.resize(CLData->number_of_clusters);
    std::cout << "Entering loop" << std::endl;
    for (int i = 0; i < CLData->number_of_clusters; i++)
    {
        centroidPoints[i] = new PointStruct;
        centroidPoints[i]->id = "";
        centroidPoints[i]->coords.resize(CLData->dimension);
    }
    std::cout << "Finished calculating centroid points..." << std::endl;

    std::vector<PointPtr> tempCentroidPoints = k_means(inputPoints, CLData->number_of_clusters, CLData->dimension);
    // Translate actual points that k_means returned to virtual centroid points
    for (int i = 0; i < CLData->number_of_clusters; i++)
    {
        for (int j = 0; j < CLData->dimension; j++)
            centroidPoints[i]->coords[j] = tempCentroidPoints[i]->coords[j];
    }
    tempCentroidPoints.clear();

    std::vector<Cluster> clusters;
    clusters.resize(CLData->number_of_clusters);
    for (int i = 0; i < CLData->number_of_clusters; i++)
    {
        clusters[i].centroidPoint = centroidPoints[i];
        clusters[i].size = 0;
    }

    std::cout << "Assigning points to clusters..." << std::endl;

    auto cluster_start = std::chrono::high_resolution_clock::now();
    execCluster(CLData, &clusters, &inputPoints, &centroidPoints);
    auto cluster_end = std::chrono::high_resolution_clock::now();

    int tCluster = std::chrono::duration_cast<std::chrono::milliseconds>(cluster_end - cluster_start).count();

    std::cout << "Evaluating silhouette..." << std::endl;
    double totalSilhouette = evalSilhouette(CLData, &clusters);

    std::cout << "Writing output file..." << std::endl;
    if (!writeToOutput(CLData, &clusters, &centroidPoints, totalSilhouette, tCluster))
    {
        std::cerr << "Error in writing output" << std::endl;
        return EXIT_FAIL_OUTPUT_ERR;
    }

    //Deleting Data Structures

    return EXIT_SUCCESS;
}