#include <iostream>
#include <algorithm>
#include <fstream>

#include "lshUtils.h"

inputData *getInputData(int *argc, char **argv)
{

    inputData *SearchData = new inputData;
    std::vector<std::string> found;
    std::string algorithmName;
    SearchData->distance_true_visible = false;

    for (int i = 0; i < *argc; i++)
    {

        if (std::string(argv[i]) == "-i")
        {
            SearchData->inputFileName = std::string(argv[i + 1]);
            std::cout << SearchData->inputFileName << std::endl;
            found.push_back("inputFile");
        }
        else if (std::string(argv[i]) == "-q")
        {
            SearchData->queryFileName = std::string(argv[i + 1]);
            std::cout << SearchData->queryFileName << std::endl;
            found.push_back("queryFile");
        }
        else if (std::string(argv[i]) == "-o")
        {
            SearchData->outputFileName = std::string(argv[i + 1]);
            std::cout << SearchData->outputFileName << std::endl;
            found.push_back("outputFile");
        }
        else if (std::string(argv[i]) == "-k")
        {
            SearchData->numberOfHyperplanes = atoi(argv[i + 1]);
            std::cout << SearchData->numberOfHyperplanes << std::endl;
            found.push_back("k");
        }
        else if (std::string(argv[i]) == "-L")
        {
            SearchData->intL = atoi(argv[i + 1]);
            std::cout << SearchData->intL << std::endl;
            found.push_back("l");
        }
        else if (std::string(argv[i]) == "-M")
        {
            SearchData->intM = atoi(argv[i + 1]);
            std::cout << SearchData->intM << std::endl;
            found.push_back("m");
        }
        else if (std::string(argv[i]) == "-probes")
        {
            SearchData->probes = atoi(argv[i + 1]);
            std::cout << SearchData->probes << std::endl;
            found.push_back("probes");
        }
        else if (std::string(argv[i]) == "-algorithm")
        {
            algorithmName = atoi(argv[i + 1]);
            std::cout << SearchData->algorithm << std::endl;
            found.push_back("algorithm");
        }
        else if (std::string(argv[i]) == "--dist-true=visible")
        {
            SearchData->distance_true_visible = true;
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
        SearchData->inputFileName = word;
    }
    if (std::find(found.begin(), found.end(), "queryFile") == found.end()) // if not found queryFile
    {
        std::cout << "Please give query file name" << std::endl;

        while ((ch = getchar()) != '\n')
            word += ch;
        SearchData->queryFileName = word;
    }
    if (std::find(found.begin(), found.end(), "outputFile") == found.end()) // if not found outputFile
    {
        std::cout << "Please give output file name" << std::endl;

        while ((ch = getchar()) != '\n')
            word += ch;
        SearchData->outputFileName = word;
    }
    if (std::find(found.begin(), found.end(), "algorithm") == found.end()) // if not found algorithm
    {
        std::cout << "Please give algorithm to be used" << std::endl;

        while ((ch = getchar()) != '\n')
            word += ch;
        algorithmName = word;
    }

    if (SearchData->inputFileName.empty() || SearchData->outputFileName.empty() || SearchData->queryFileName.empty() || algorithmName.empty())
    {
        std::cerr << "The name of the input file, output file, query file and algorithm must be given."
                  << std::endl;
        delete SearchData;
        return NULL;
    }

    if (algorithmName == "LSH")
        SearchData->algorithm = ALG_LSH;
    else if (algorithmName == "Hypercube")
        SearchData->algorithm = ALG_HC;
    else if (algorithmName == "Frechet")
        SearchData->algorithm = ALG_FR;
    else
    {
        std::cerr << "Unknown algorithm. Valid options are LSH, Hypercube, or Frechet" << std::endl;
        delete SearchData;
        return NULL;
    }

    return SearchData;
}

int writeToOutput(inputData *LSHData,
                  std::vector<PointPtr> queryPoints,
                  std::vector<kNeighboursPtr> queryOutputData,
                  std::vector<kNeighboursPtr> queryTrueNeighbors,
                  std::vector<std::vector<PointPtr>> queryRangeSearch,
                  std::vector<double> tLSH,
                  std::vector<double> tTrue)
{
    std::ofstream outputFile(LSHData->outputFileName);
    if (!outputFile.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << LSHData->outputFileName << "'"
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
                       << "distanceLSH: " << queryOutputData[i]->neighbours[j]->dist << std::endl;
            if (LSHData->distance_true_visible)
            {
                outputFile << "True Nearest neighbor-"
                           << j + 1 << ": " << queryTrueNeighbors[i]->neighbours[j]->point->id << std::endl;
            }
            outputFile << "distanceTrue: " << queryTrueNeighbors[i]->neighbours[j]->dist << std::endl;
        }

        outputFile << "tLSH: " << (double)(tLSH[i] / 1000) << 's' << std::endl
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
                inputData *LSHData)
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
    delete LSHData;
}