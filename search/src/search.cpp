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
#include "lsh_frechet_dsc.h"
#include "lsh_frechet_cont.h"
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
        inputData *SearchData = getInputData(&argc, argv);
        if (SearchData == NULL)
            return EXIT_FAILURE;

        // Opening inputFile
        std::cout << "Reading input file " << SearchData->inputFileName << "..." << std::endl;
        std::vector<std::string> inputLines = get_lines(SearchData->inputFileName);

        if (SearchData->algorithm == ALG_LSH)
        {

            // Getting points from lines
            std::vector<PointPtr> inputPoints;
            SearchData->dimension = get_points(inputLines, &inputPoints);
            int numOfInputPoints = inputPoints.size();

            HashTables HashTablesObject(SearchData->intL, SearchData->numberOfHyperplanes, numOfInputPoints, SearchData->dimension, numOfInputPoints / 4);
            // Inserting points to hash table
            std::cout << "Inserting items to hash table..." << std::endl;
            for (int i = 0; i < inputPoints.size(); i++)
                HashTablesObject.HashTables::InsertPoint(inputPoints[i]);

            // Getting lines from query file
            std::cout << "Reading query file " << SearchData->queryFileName << "..." << std::endl;
            std::vector<std::string> queryLines = get_lines(SearchData->queryFileName);

            // Getting points from lines
            std::vector<PointPtr> queryPoints;
            get_points(queryLines, &queryPoints);

            // LSH k nearest neighbor search
            std::vector<std::vector<Neighbour> *> k_nearest_neighbours;
            int numOfQueries = queryLines.size();
            k_nearest_neighbours.resize(numOfQueries);

            std::vector<kNeighboursPtr> queryOutputData;
            queryOutputData.resize(numOfQueries);

            std::vector<double> tLSH;
            tLSH.resize(numOfQueries);

            std::cout << "Executing LSH search algorithm..." << std::endl;

            for (int i = 0; i < numOfQueries; i++)
            {
                auto LSH_start = std::chrono::high_resolution_clock::now();
                queryOutputData[i] = HashTablesObject.HashTables::find_k_nearest_neighbours(queryPoints[i], 1);
                auto LSH_end = std::chrono::high_resolution_clock::now();
                tLSH[i] = std::chrono::duration_cast<std::chrono::milliseconds>(LSH_end - LSH_start).count();
            }
            double tLSHAverage = 0.0;
            for (int i = 0; i < numOfQueries; i++)
                tLSHAverage += tLSH[i];

            // Brute force k nearest neighbor search
            std::vector<kNeighboursPtr>
                queryTrueNeighbors;
            queryTrueNeighbors.resize(numOfQueries);

            std::vector<double> tTrue;
            tTrue.resize(numOfQueries);

            std::cout << "Executing brute-force search algorithm..." << std::endl;

            for (int i = 0; i < numOfQueries; i++)
            {
                auto True_start = std::chrono::high_resolution_clock::now();
                queryTrueNeighbors[i] = find_k_true_neighbours(queryPoints[i], 1, inputPoints, SearchData->dimension);
                auto True_end = std::chrono::high_resolution_clock::now();
                tTrue[i] = std::chrono::duration_cast<std::chrono::milliseconds>(True_end - True_start).count();
            }
            double tTrueAverage = 0.0;
            for (int i = 0; i < numOfQueries; i++)
                tTrueAverage += tLSH[i];

            // Writing results to outputFile
            if (writeToOutput(SearchData, queryPoints, queryOutputData, queryTrueNeighbors, tLSHAverage, tTrueAverage, "LSH_Vector"))
                return EXIT_FAILURE;

            // Deleting Data Structures
            std::cout << "Freeing memory" << std::endl;
            deleteData(&inputPoints, &queryPoints, &k_nearest_neighbours, &queryOutputData, &queryTrueNeighbors, SearchData);
        }
        else if (SearchData->algorithm == ALG_FR)
        {
            if (SearchData->metric == MTR_DISC)
            {

                // Getting points from lines
                std::vector<PointPtr> inputPoints;
                SearchData->dimension = get_points(inputLines, &inputPoints);
                int numOfInputPoints = inputPoints.size();
                std::cout << "Dimension: " << SearchData->dimension << std::endl;

                // std::vector<PointPtr> snappedPoints;
                // for (int i = 0; i < numOfInputPoints; i++)
                //     snappedPoints.push_back(snap_point(inputPoints[i], DELTA, SearchData->dimension));

                // std::vector<PointPtr> removedDupPoints;

                // // Turning points from lines into 2d-projected curves
                // std::vector<CurvePtr> *curves;
                // SearchData->dimension = get_curves(inputLines, curves);
                // int numOfInputPoints = curves->size();

                // Create HashTablesObject that stores curves, projected to vectors of size 2*dim
                FrechetDiscreteHashTables HashTablesObject(SearchData->intL, SearchData->numberOfHyperplanes, numOfInputPoints, SearchData->dimension, numOfInputPoints / 4);
                // Inserting curves to hash table, after snapping them
                std::cout << "Inserting items to hash table..." << std::endl;
                for (int i = 0; i < numOfInputPoints; i++)
                {
                    std::cout << inputPoints[i]->coords.size() << std::endl;
                    HashTablesObject.FrechetDiscreteHashTables::FrDscInsertPoint(inputPoints[i]);
                }
                // Getting lines from query file
                std::cout << "Reading query file " << SearchData->queryFileName << "..." << std::endl;
                std::vector<std::string> queryLines = get_lines(SearchData->queryFileName);
                // Getting curves from query lines
                std::vector<PointPtr> queryCurves;
                std::cout << "get_curves got " << get_curves(queryLines, &queryCurves) << std::endl;
                int numOfQueries = queryLines.size();
                // LSH k nearest neighbor search
                std::vector<std::vector<Neighbour> *> k_nearest_neighbours;
                k_nearest_neighbours.resize(numOfQueries);

                kNeighboursPtr knnOutputData;

                std::vector<NeighbourPtr> queryOutputData;
                queryOutputData.resize(numOfQueries);

                std::vector<double> tLSH;
                tLSH.resize(numOfQueries);

                std::cout << "Executing LSH search algorithm..." << std::endl;

                for (int i = 0; i < numOfQueries; i++)
                {
                    auto LSH_start = std::chrono::high_resolution_clock::now();
                    knnOutputData = HashTablesObject.FrechetDiscreteHashTables::FrDsc_find_k_nearest_neighbours(queryCurves[i], 1);
                    queryOutputData[i] = knnOutputData->neighbours[0];
                    auto LSH_end = std::chrono::high_resolution_clock::now();
                    tLSH[i] = std::chrono::duration_cast<std::chrono::milliseconds>(LSH_end - LSH_start).count();
                }
                double tLSHAverage = 0.0;
                for (int i = 0; i < numOfQueries; i++)
                    tLSHAverage += tLSH[i];

                // Brute force k nearest neighbor search

                kNeighboursPtr kTrueOutputData;
                std::vector<NeighbourPtr> queryTrueNeighbors;
                queryTrueNeighbors.resize(numOfQueries);

                std::vector<double> tTrue;
                tTrue.resize(numOfQueries);

                std::cout << "Executing brute-force search algorithm..." << std::endl;

                for (int i = 0; i < numOfQueries; i++)
                {
                    auto True_start = std::chrono::high_resolution_clock::now();
                    kTrueOutputData = find_k_true_neighbours(queryCurves[i], 1, inputPoints, SearchData->dimension);
                    queryTrueNeighbors[i] = kTrueOutputData->neighbours[0];
                    auto True_end = std::chrono::high_resolution_clock::now();
                    tTrue[i] = std::chrono::duration_cast<std::chrono::milliseconds>(True_end - True_start).count();
                }
                double tTrueAverage = 0.0;
                for (int i = 0; i < numOfQueries; i++)
                    tTrueAverage += tLSH[i];
                
                // LSH range search
                std::vector<std::vector<PointPtr>> queryRangeSearch;
                queryRangeSearch.resize(queryLines.size());
                std::cout << "Executing range search algorithm..." << std::endl;
                for (int i = 0; i < queryLines.size(); i++)
                {
                    queryRangeSearch[i] = HashTablesObject.FrechetDiscreteHashTables::range_search(queryCurves[i], 100);
                }

                // Writing results to outputFile
                if (writeToOutputFrDsc(SearchData, queryCurves, queryOutputData, queryTrueNeighbors, tLSHAverage, tTrueAverage, "LSH_Vector"))
                    return EXIT_FAILURE;

                // Deleting Data Structures
                // std::cout << "Freeing memory" << std::endl;
                // deleteData(&inputPoints, &queryPoints, &k_nearest_neighbours, &queryOutputData, &queryTrueNeighbors, SearchData);
            }
            else if (SearchData->metric == MTR_CONT)
            {
                // Getting points from lines
                std::vector<PointPtr> inputPoints;
                SearchData->dimension = get_points(inputLines, &inputPoints);
                int numOfInputPoints = inputPoints.size();

                // Create HashTablesObject that stores curves, projected to vectors of size 2*dim
                FrechetContinuousHashTables HashTablesObject(1, SearchData->numberOfHyperplanes, numOfInputPoints, SearchData->dimension, numOfInputPoints / 4);

                // Inserting curves to hash table, after snapping them
                std::cout << "Inserting items to hash table..." << std::endl;
                for (int i = 0; i < numOfInputPoints; i++)
                    HashTablesObject.FrechetContinuousHashTables::FrContInsertPoint(inputPoints[i]);

                // Getting lines from query file
                std::cout << "Reading query file " << SearchData->queryFileName << "..." << std::endl;
                std::vector<std::string> queryLines = get_lines(SearchData->queryFileName);

                // Getting curves from query lines
                std::vector<PointPtr> queryCurves;
                get_curves(queryLines, &queryCurves);

                // LSH k nearest neighbor search
                std::vector<std::vector<Neighbour> *> k_nearest_neighbours;
                int numOfQueries = queryLines.size();
                k_nearest_neighbours.resize(numOfQueries);

                kNeighboursPtr knnOutputData;

                std::vector<NeighbourPtr> queryOutputData;
                queryOutputData.resize(numOfQueries);

                std::vector<double> tLSH;
                tLSH.resize(numOfQueries);

                std::cout << "Executing LSH search algorithm..." << std::endl;

                for (int i = 0; i < numOfQueries; i++)
                {
                    auto LSH_start = std::chrono::high_resolution_clock::now();
                    knnOutputData = HashTablesObject.FrechetContinuousHashTables::FrCont_find_k_nearest_neighbours(queryCurves[i], 1);
                    queryOutputData[i] = knnOutputData->neighbours[0];
                    auto LSH_end = std::chrono::high_resolution_clock::now();
                    tLSH[i] = std::chrono::duration_cast<std::chrono::milliseconds>(LSH_end - LSH_start).count();
                }
                double tLSHAverage = 0.0;
                for (int i = 0; i < numOfQueries; i++)
                    tLSHAverage += tLSH[i];

                // Brute force k nearest neighbor search

                kNeighboursPtr kTrueOutputData;
                std::vector<NeighbourPtr> queryTrueNeighbors;
                queryTrueNeighbors.resize(numOfQueries);

                std::vector<double> tTrue;
                tTrue.resize(numOfQueries);

                std::cout << "Executing brute-force search algorithm..." << std::endl;

                for (int i = 0; i < numOfQueries; i++)
                {
                    auto True_start = std::chrono::high_resolution_clock::now();
                    kTrueOutputData = find_k_true_neighbours(queryCurves[i], 1, inputPoints, SearchData->dimension);
                    queryTrueNeighbors[i] = kTrueOutputData->neighbours[0];
                    auto True_end = std::chrono::high_resolution_clock::now();
                    tTrue[i] = std::chrono::duration_cast<std::chrono::milliseconds>(True_end - True_start).count();
                }
                double tTrueAverage = 0.0;
                for (int i = 0; i < numOfQueries; i++)
                    tTrueAverage += tLSH[i];

                // Writing results to outputFile
                if (writeToOutputFrDsc(SearchData, queryCurves, queryOutputData, queryTrueNeighbors, tLSHAverage, tTrueAverage, "LSH_Vector"))
                    return EXIT_FAILURE;

                // Deleting Data Structures
                // std::cout << "Freeing memory" << std::endl;
                // deleteData(&inputPoints, &queryPoints, &k_nearest_neighbours, &queryOutputData, &queryTrueNeighbors, SearchData);
            }
        }
        return EXIT_SUCCESS;

        option = checkRerun();
    }
    return EXIT_SUCCESS;
}