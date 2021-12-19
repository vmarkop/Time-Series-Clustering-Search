#include <iostream>
#include <algorithm>
#include <fstream>

#include "lshUtils.h"

#define DEF_M 10
#define DEF_PROBES 2

inputData *getInputData(int *argc, char **argv)
{

    inputData *SearchData = new inputData;
    std::vector<std::string> found;
    std::string algorithmName, metricName;

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
            algorithmName = std::string(argv[i + 1]);
            std::cout << algorithmName << std::endl;
            found.push_back("algorithm");
        }
        else if (std::string(argv[i]) == "-metric")
        {
            metricName = std::string(argv[i + 1]);
            std::cout << metricName << std::endl;
            found.push_back("metric");
        }
    }

    //-----------------------------------------------------------------------------------------//

    char ch;
    std::string word = "";

    *argc = 1;

    found.push_back(" ");
    std::string input = {};

    if (std::find(found.begin(), found.end(), "k") == found.end()) // if not found -k
    {
        SearchData->numberOfHyperplanes = DEF_K;
    }
    if (std::find(found.begin(), found.end(), "l") == found.end()) // if not found -l
    {
        SearchData->intL = DEF_L;
    }
    if (std::find(found.begin(), found.end(), "m") == found.end()) // if not found -m
    {
        SearchData->intM = DEF_M;
    }
    if (std::find(found.begin(), found.end(), "probes") == found.end()) // if not found -probes
    {
        SearchData->intM = DEF_PROBES;
    }

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
        word = "";
        while ((ch = getchar()) != '\n')
            word += ch;
        SearchData->queryFileName = word;
    }
    if (std::find(found.begin(), found.end(), "outputFile") == found.end()) // if not found outputFile
    {
        std::cout << "Please give output file name" << std::endl;
        word = "";
        while ((ch = getchar()) != '\n')
            word += ch;
        SearchData->outputFileName = word;
    }
    if (std::find(found.begin(), found.end(), "algorithm") == found.end()) // if not found algorithm
    {
        std::cout << "Please give algorithm to be used" << std::endl;
        word = "";
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
    {
        SearchData->algorithm = ALG_FR;
        if (std::find(found.begin(), found.end(), "metric") == found.end()) // if not found metric
        {
            std::cout << "Please give metric: discrete or continuous to be used" << std::endl;
            word = "";
            while ((ch = getchar()) != '\n')
                word += ch;
            metricName = word;
        }

        if (metricName == "discrete")
            SearchData->metric = MTR_DISC;
        else if (metricName == "continuous")
            SearchData->metric = MTR_CONT;
        else
        {
            std::cerr << "Unknown Frechet metric. Valid options are discrete or continuous." << std::endl;
            delete SearchData;
            return NULL;
        }
    }
    else
    {
        std::cerr << "Unknown algorithm. Valid options are LSH, Hypercube, or Frechet" << std::endl;
        delete SearchData;
        return NULL;
    }

    return SearchData;
}

int writeToOutput(inputData *SearchData,
                  std::vector<PointPtr> queryPoints,
                  std::vector<kNeighboursPtr> queryOutputData,
                  std::vector<kNeighboursPtr> queryTrueNeighbors,
                  double tLSH,
                  double tTrue,
                  std::string algorithm)
{
    std::cout << "Creating output file..." << std::endl;
    std::ofstream outputFile(SearchData->outputFileName);
    if (!outputFile.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << SearchData->outputFileName << "'"
                  << std::endl;
        return EXIT_FAILURE;
    }

    for (int i = 0; i < queryPoints.size(); i++)
    {
        outputFile << "Query: "
                   << queryPoints[i]->id << std::endl
                   << "Algorithm: " << algorithm << std::endl;

        outputFile << "Approximate Nearest neighbor: "
                   << queryOutputData[i]->neighbours[0]->point->id << std::endl;

        outputFile << "True Nearest neighbor: "
                   << queryTrueNeighbors[i]->neighbours[0]->point->id << std::endl;

        outputFile << "distanceApproximate: "
                   << queryOutputData[i]->neighbours[0]->dist << std::endl;

        outputFile << "distanceTrue: "
                   << queryTrueNeighbors[i]->neighbours[0]->dist << std::endl;

        outputFile << std::endl;
    }

    outputFile << "tApproximateAverage: " << (double)(tLSH / queryPoints.size() / 1000) << 's' << std::endl
               << "tTrueAverage: " << (double)(tTrue / queryPoints.size() / 1000) << 's' << std::endl;

    outputFile.close();
    return EXIT_SUCCESS;
}

int writeToOutputFrDsc(inputData *SearchData,
                       std::vector<PointPtr> queryPoints,
                       std::vector<NeighbourPtr> queryOutputData,
                       std::vector<NeighbourPtr> queryTrueNeighbors,
                       double tLSH,
                       double tTrue,
                       std::string algorithm)
{
    std::cout << "Creating output file..." << std::endl;
    std::ofstream outputFile(SearchData->outputFileName);
    if (!outputFile.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << SearchData->outputFileName << "'"
                  << std::endl;
        return EXIT_FAILURE;
    }

    for (int i = 0; i < queryPoints.size(); i++)
    {
        outputFile << "Query: "
                   << queryPoints[i]->id << std::endl
                   << "Algorithm: " << algorithm << std::endl;

        outputFile << "Approximate Nearest neighbor: "
                   << queryOutputData[i]->point->id << std::endl;

        outputFile << "True Nearest neighbor: "
                   << queryTrueNeighbors[i]->point->id << std::endl;

        outputFile << "distanceApproximate: "
                   << queryOutputData[i]->dist << std::endl;

        outputFile << "distanceTrue: "
                   << queryTrueNeighbors[i]->dist << std::endl;

        outputFile << std::endl;
    }

    outputFile << "tApproximateAverage: " << (double)(tLSH / queryPoints.size() / 1000) << 's' << std::endl
               << "tTrueAverage: " << (double)(tTrue / queryPoints.size() / 1000) << 's' << std::endl;

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