#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <climits>

#include "hashTable.h"
#include "mathUtils.h"
#include "lshUtils.h"

int main(int argc, char **argv)
{
    std::string option;
    while (!(option == "TERM" || option == "term"))
    {
        std::cout << std::endl
                  << std::endl
                  << "STARTING" << std::endl;

        // Reading arguments
        inputData *LSHData = getInputData(&argc, argv);
        if (LSHData == NULL)
            return EXIT_FAILURE;

        // Opening inputFile
        std::cout << "Reading input file " << LSHData->inputFileName << "..." << std::endl;
        std::vector<std::string> inputLines = get_lines(LSHData->inputFileName);

        // Getting points from lines
        std::vector<PointPtr> inputPoints;
        LSHData->dimension = get_points(inputLines, &inputPoints);

        std::cout << "Dimension: " << LSHData->dimension << std::endl
                  << std::endl;

        int numOfInputPoints = inputPoints.size();
        HashTables HashTablesObject(LSHData->intL, LSHData->numberOfHyperplanes, numOfInputPoints, LSHData->dimension, numOfInputPoints / 4);

        // Inserting points to hash table
        std::cout << "Inserting items to hash table..." << std::endl;
        for (int i = 0; i < inputPoints.size(); i++)
            HashTablesObject.HashTables::InsertPoint(inputPoints[i]);

        // Getting lines from query file
        std::cout << "Reading query file " << LSHData->queryFileName << "..." << std::endl;
        std::vector<std::string> queryLines = get_lines(LSHData->queryFileName);

        // Getting points from lines
        std::vector<PointPtr> queryPoints;
        get_points(queryLines, &queryPoints);

        // LSH k nearest neighbor search
        std::vector<std::vector<Neighbour> *> k_nearest_neighbours;
        k_nearest_neighbours.resize(queryLines.size());

        std::vector<kNeighboursPtr> queryOutputData;
        queryOutputData.resize(queryLines.size());

        std::vector<double> tLSH;
        tLSH.resize(queryLines.size());

        std::cout << "Executing LSH search algorithm..." << std::endl;

        for (int i = 0; i < queryLines.size(); i++)
        {

            auto LSH_start = std::chrono::high_resolution_clock::now();
            queryOutputData[i] = HashTablesObject.HashTables::find_k_nearest_neighbours(queryPoints[i], LSHData->numberOfNearest);
            auto LSH_end = std::chrono::high_resolution_clock::now();
            tLSH[i] = std::chrono::duration_cast<std::chrono::milliseconds>(LSH_end - LSH_start).count();
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
            queryTrueNeighbors[i] = find_k_true_neighbours(queryPoints[i], LSHData->numberOfNearest, inputPoints, LSHData->dimension);
            auto True_end = std::chrono::high_resolution_clock::now();
            tTrue[i] = std::chrono::duration_cast<std::chrono::milliseconds>(True_end - True_start).count();
        }

        // LSH range search
        std::vector<std::vector<PointPtr>> queryRangeSearch;
        queryRangeSearch.resize(queryLines.size());

        std::cout << "Executing range search algorithm..." << std::endl;

        for (int i = 0; i < queryLines.size(); i++)
        {
            queryRangeSearch[i] = HashTablesObject.HashTables::range_search(queryPoints[i], LSHData->radius);
        }

        // Writing results to outputFile
        if (writeToOutput(LSHData, queryPoints, queryOutputData, queryTrueNeighbors, queryRangeSearch, tLSH, tTrue))
            return EXIT_FAILURE;

        // Deleting Data Structures
        deleteData(&inputPoints, &queryPoints, &k_nearest_neighbours, &queryOutputData, &queryTrueNeighbors, LSHData);

        option = checkRerun();
    }
    return EXIT_SUCCESS;
}