#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <climits>

#include "hypercubeUtils.h"
#include "HChashTable.h"
#include "mathUtils.h"

int main(int argc, char **argv)
{
    std::string option;
    while (!(option == "TERM" || option == "term"))
    {
        std::cout << std::endl
                  << std::endl
                  << "STARTING" << std::endl;

        // Reading arguments
        inputData *HCData = getInputData(&argc, argv);
        if (HCData == NULL)
            return EXIT_FAILURE;

        // Opening inputFile
        std::cout << "Reading input file " << HCData->inputFileName << "..." << std::endl;
        std::vector<std::string> inputLines = get_lines(HCData->inputFileName);

        // Getting points from lines
        std::vector<PointPtr> inputPoints;
        HCData->dimension = get_points(inputLines, &inputPoints);

        std::cout << "Dimension:" << HCData->dimension << std::endl
                  << std::endl;

        HChashTable HypercubeObject(HCData->dimension, HCData->projectionDimension, HCData->probes, HCData->maxCandidatePoints);

        // Inserting points to hash table
        std::cout << "Inserting items to hash table..." << std::endl;
        for (int i = 0; i < inputPoints.size(); i++)
            HypercubeObject.HChashTable::InsertPoint(inputPoints[i]);

        // Getting lines from query file
        std::cout << "Reading query file " << HCData->queryFileName << "..." << std::endl;
        std::vector<std::string> queryLines = get_lines(HCData->queryFileName);

        // Getting points from lines
        std::vector<PointPtr> queryPoints;
        get_points(queryLines, &queryPoints);

        // HyperCube k nearest neighbor search
        std::vector<std::vector<Neighbour> *> k_nearest_neighbours;
        k_nearest_neighbours.resize(queryLines.size());

        std::vector<kNeighboursPtr> queryOutputData;
        queryOutputData.resize(queryLines.size());

        std::vector<double> tCube;
        tCube.resize(queryLines.size());

        std::cout << "Executing Hypercube search algorithm..." << std::endl;

        for (int i = 0; i < queryLines.size(); i++)
        {
            auto Cube_start = std::chrono::high_resolution_clock::now();
            queryOutputData[i] = HypercubeObject.HChashTable::find_k_nearest_neighbours(queryPoints[i], HCData->numberOfNearest);
            auto Cube_end = std::chrono::high_resolution_clock::now();
            tCube[i] = std::chrono::duration_cast<std::chrono::milliseconds>(Cube_end - Cube_start).count();
        }

        // Brute force k nearest neighbor search
        std::vector<kNeighboursPtr> queryTrueNeighbors;
        queryTrueNeighbors.resize(queryLines.size());

        std::vector<double> tTrue;
        tTrue.resize(queryLines.size());

        std::cout << "Executing brute-force search algorithm..." << std::endl;

        for (int i = 0; i < queryLines.size(); i++)
        {
            auto True_start = std::chrono::high_resolution_clock::now();
            queryTrueNeighbors[i] = find_k_true_neighbours(queryPoints[i], HCData->numberOfNearest, inputPoints, HCData->dimension);
            auto True_end = std::chrono::high_resolution_clock::now();
            tTrue[i] = std::chrono::duration_cast<std::chrono::milliseconds>(True_end - True_start).count();
        }

        // HyperCube range search
        std::vector<std::vector<PointPtr>> queryRangeSearch;
        queryRangeSearch.resize(queryLines.size());

        std::cout << "Executing range search algorithm..." << std::endl;

        for (int i = 0; i < queryLines.size(); i++)
        {
            queryRangeSearch[i] = HypercubeObject.HChashTable::range_search(queryPoints[i], HCData->radius);
        }

        // Writing results to outputFile
        if (writeToOutput(HCData, queryPoints, queryOutputData, queryTrueNeighbors, queryRangeSearch, tCube, tTrue))
            return EXIT_FAILURE;

        // Deleting Data Structures
        deleteData(&inputPoints, &queryPoints, &k_nearest_neighbours, &queryOutputData, &queryTrueNeighbors, HCData);

        option = checkRerun();
    }
    return EXIT_SUCCESS;
}