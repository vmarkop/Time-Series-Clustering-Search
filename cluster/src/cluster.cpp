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
    std::cout << "Initialization..." << std::endl;
    std::vector<PointPtr> centroidPoints;
    centroidPoints.resize(CLData->number_of_clusters);
    for (int i = 0; i < CLData->number_of_clusters; i++)
    {
        centroidPoints[i] = new PointStruct;
        centroidPoints[i]->id = "";
        centroidPoints[i]->coords.resize(CLData->dimension);
    }

    std::vector<PointPtr> tempCentroidPoints = k_means(inputPoints, CLData); // CLData->number_of_clusters, CLData->dimension);
    // Translate actual points that k_means returned to virtual centroid points
    for (int i = 0; i < CLData->number_of_clusters; i++)
    {
        if (CLData->update == UPDATE_FRECHET) // turn centroid points into 2d points
        {
            PointPtr concated_point = concat_point(tempCentroidPoints[i], CLData->dimension);
            centroidPoints[i] = concated_point;
        }
        else
        {
            for (int j = 0; j < CLData->dimension; j++)
                centroidPoints[i]->coords[j] = tempCentroidPoints[i]->coords[j];
        }
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

    std::vector<PointPtr> *inputPoints_2d = new std::vector<PointPtr>;
    inputPoints_2d->resize(CLData->numberOfInputPoints);
    for (int i = 0; i < CLData->numberOfInputPoints; i++)
    {
        (*inputPoints_2d)[i] = concat_point(inputPoints[i], CLData->dimension);
    }

    auto cluster_start = std::chrono::high_resolution_clock::now();
    std::cout << "Assigning points to clusters...done" << std::endl;
    if (execCluster(CLData, &clusters, &inputPoints, inputPoints_2d, &centroidPoints) == EXIT_FAILURE)
    {
        return EXIT_FAILURE;
    }
    auto cluster_end = std::chrono::high_resolution_clock::now();

    int tCluster = std::chrono::duration_cast<std::chrono::milliseconds>(cluster_end - cluster_start).count();

    std::cout << "Evaluating silhouette..." << std::endl;
    double totalSilhouette = evalSilhouette(CLData, &clusters);

    std::cout << "Writing output file..." << std::endl;
    if (writeToOutput(CLData, &clusters, &centroidPoints, totalSilhouette, tCluster))
    {
        std::cerr << "Error in writing output" << std::endl;
        return EXIT_FAIL_OUTPUT_ERR;
    }

    // Deleting Data Structures

    return EXIT_SUCCESS;
}