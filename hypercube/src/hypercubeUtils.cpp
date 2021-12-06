#include <iostream>
#include <time.h>
#include <algorithm>
#include <stdlib.h>
#include <fstream>
#include <random>

#include "hypercubeUtils.h"

inputData *getInputData(int *argc, char **argv)
{

    inputData *HCData = new inputData;
    std::vector<std::string> found;

    for (int i = 0; i < *argc; i++)
    {

        if (std::string(argv[i]) == "-i")
        {
            HCData->inputFileName = std::string(argv[i + 1]);
            std::cout << HCData->inputFileName << std::endl;
            found.push_back("inputFile");
        }
        else if (std::string(argv[i]) == "-q")
        {
            HCData->queryFileName = std::string(argv[i + 1]);
            std::cout << HCData->queryFileName << std::endl;
            found.push_back("queryFile");
        }
        else if (std::string(argv[i]) == "-o")
        {
            HCData->outputFileName = std::string(argv[i + 1]);
            std::cout << HCData->outputFileName << std::endl;
            found.push_back("outputFile");
        }
        else if (std::string(argv[i]) == "-k")
        {
            HCData->projectionDimension = atoi(argv[i + 1]);
            std::cout << HCData->projectionDimension << std::endl;
            found.push_back("k");
        }
        else if (std::string(argv[i]) == "-M")
        {
            HCData->maxCandidatePoints = atoi(argv[i + 1]);
            std::cout << HCData->maxCandidatePoints << std::endl;
            found.push_back("m");
        }
        else if (std::string(argv[i]) == "-N")
        {
            HCData->numberOfNearest = atoi(argv[i + 1]);
            std::cout << HCData->numberOfNearest << std::endl;
            found.push_back("n");
        }
        else if (std::string(argv[i]) == "-probes")
        {
            HCData->probes = atoi(argv[i + 1]);
            std::cout << HCData->probes << std::endl;
            found.push_back("probes");
        }
        else if (std::string(argv[i]) == "-R")
        {
            HCData->radius = atoi(argv[i + 1]);
            std::cout << HCData->radius << std::endl;
            found.push_back("r");
        }
        else if (std::string(argv[i]) == "--dist-true=visible")
        {
        }
    }

    //-----------------------------------------------------------------------------------------//

    char ch;
    std::string word = "";

    *argc = 1;

    found.push_back(" ");
    std::string input = {};

    if (std::find(found.begin(), found.end(), "inputFile") == found.end()) // if not found inputFile
    {
        std::cout << "Please give input file name" << std::endl;

        while ((ch = getchar()) != '\n')
            word += ch;
        HCData->inputFileName = word;
    }
    if (std::find(found.begin(), found.end(), "queryFile") == found.end()) // if not found queryFile
    {
        std::cout << "Please give query file name" << std::endl;

        while ((ch = getchar()) != '\n')
            word += ch;
        HCData->queryFileName = word;
    }
    if (std::find(found.begin(), found.end(), "outputFile") == found.end()) // if not found outputFile
    {
        std::cout << "Please give output file name" << std::endl;

        while ((ch = getchar()) != '\n')
            word += ch;
        HCData->outputFileName = word;
    }
    if (std::find(found.begin(), found.end(), "k") == found.end()) // if not found inputFile
    {
        std::cout << "Please give k: Press [ENTER] for default value or type the desired value." << std::endl;
        HCData->projectionDimension = DEF_K;
        word = "";
        while ((ch = getchar()) != '\n')
            word += ch;
        if (word == "")
            std::cout << "Using Default Value of k = " << HCData->projectionDimension << std::endl;
        else
            HCData->projectionDimension = stoi(word);
    }
    if (std::find(found.begin(), found.end(), "m") == found.end()) // if not found inputFile
    {
        std::cout << "Please give M: Press [ENTER] for default value or type the desired value." << std::endl;
        HCData->maxCandidatePoints = DEF_M;
        word = "";
        while ((ch = getchar()) != '\n')
            word += ch;
        if (word == "")
            std::cout << "Using Default Value of M = " << HCData->maxCandidatePoints << std::endl;
        else
            HCData->maxCandidatePoints = stoi(word);
    }
    if (std::find(found.begin(), found.end(), "probes") == found.end()) // if not found inputFile
    {
        std::cout << "Please give number of probes: Press [ENTER] for default value or type the desired value." << std::endl;
        HCData->probes = DEF_PROBES;
        word = "";
        while ((ch = getchar()) != '\n')
            word += ch;
        if (word == "")
            std::cout << "Using Default Value of probes = " << HCData->probes << std::endl;
        else
            HCData->probes = stoi(word);
    }
    if (std::find(found.begin(), found.end(), "n") == found.end()) // if not found inputFile
    {
        std::cout << "Please give N: Press [ENTER] for default value or type the desired value." << std::endl;
        HCData->numberOfNearest = DEF_N;
        word = "";
        while ((ch = getchar()) != '\n')
            word += ch;
        if (word == "")
            std::cout << "Using Default Value of N = " << HCData->numberOfNearest << std::endl;
        else
            HCData->numberOfNearest = stoi(word);
    }
    if (std::find(found.begin(), found.end(), "r") == found.end()) // if not found inputFile
    {
        std::cout << "Please give R: Press [ENTER] for default value or type the desired value." << std::endl;
        HCData->radius = DEF_R;
        word = "";
        while ((ch = getchar()) != '\n')
            word += ch;
        if (word == "")
            std::cout << "Using Default Value of radius = " << HCData->radius << std::endl;
        else
            HCData->radius = stoi(word);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////
    if (HCData->inputFileName.empty() || HCData->outputFileName.empty() || HCData->queryFileName.empty())
    {
        std::cerr << "Arguments must contain all input file, output file and query file. The rest of the arguments are optional"
                  << std::endl;
        return NULL;
    }
    return HCData;
}

int writeToOutput(inputData *HCData,
                  std::vector<PointPtr> queryPoints,
                  std::vector<kNeighboursPtr> queryOutputData,
                  std::vector<kNeighboursPtr> queryTrueNeighbors,
                  std::vector<std::vector<PointPtr>> queryRangeSearch,
                  std::vector<double> tCube,
                  std::vector<double> tTrue)
{
    std::ofstream outputFile(HCData->outputFileName);
    if (!outputFile.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << HCData->outputFileName << "'"
                  << std::endl;
        return EXIT_FAILURE;
    }

    for (int i = 0; i < queryPoints.size(); i++)
    {
        outputFile << "Query: "
                   << queryPoints[i]->id << std::endl;

        for (int j = 0; j < queryOutputData[i]->size; j++)
        {
            outputFile << "Nearest neighbor-"
                       << j + 1 << ": " << queryOutputData[i]->neighbours[j]->point->id << std::endl
                       << "distanceHypercube: " << queryOutputData[i]->neighbours[j]->dist << std::endl;
            // if (HCData->distance_true_visible)
            // {
            outputFile << "True Nearest neighbor-"
                       << j + 1 << ": " << queryTrueNeighbors[i]->neighbours[j]->point->id << std::endl;
            // }
            outputFile << "distanceTrue: " << queryTrueNeighbors[i]->neighbours[j]->dist << std::endl;
        }

        outputFile << "tCube: " << (double)(tCube[i] / 1000) << 's' << std::endl
                   << "tTrue: " << (double)(tTrue[i] / 1000) << 's' << std::endl
                   << "R-near neighbors:" << std::endl;
        for (int j = 0; j < queryRangeSearch[i].size(); j++)
            outputFile << queryRangeSearch[i][j]->id << std::endl;
        outputFile << std::endl
                   << std::endl;
    }
    outputFile.close();
    return EXIT_SUCCESS;
}

void deleteData(std::vector<PointPtr> *inputPoints,
                std::vector<PointPtr> *queryPoints,
                std::vector<std::vector<Neighbour> *> *k_nearest_neighbours,
                std::vector<kNeighboursPtr> *queryOutputData,
                std::vector<kNeighboursPtr> *queryTrueNeighbors,
                inputData *HCData)
{
    for (int i = 0; i < inputPoints->size(); i++)
    {
        delete (*inputPoints)[i];
        if (i < queryPoints->size()) // Query points will always be <= input points, so this is safe
        {
            delete (*queryPoints)[i];
            delete (*k_nearest_neighbours)[i];
            for (int j = 0; j < (*queryOutputData)[i]->size; j++)
            {
                delete (*queryOutputData)[i]->neighbours[j];
            }
            delete (*queryOutputData)[i];

            for (int j = 0; j < (*queryTrueNeighbors)[i]->size; j++)
            {
                delete (*queryTrueNeighbors)[i]->neighbours[j];
            }
            delete (*queryTrueNeighbors)[i];
        }
    }
    delete HCData;
}