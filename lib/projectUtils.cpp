#include <fstream>

#include "projectUtils.h"
#include "mathUtils.h"

std::vector<std::string> get_lines(std::string fileName)
{
    std::ifstream file(fileName);
    if (!file.is_open())
    {
        std::cerr << "Could not open the file: '"
                  << fileName << "'"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    std::cout << "Reading input file " << fileName << "..." << std::endl;
    std::vector<std::string> inputLines;
    std::string line;
    while (std::getline(file, line))
    {
        inputLines.push_back(line);
    }
    file.close();
    return inputLines;
}

int get_points(std::vector<std::string> linesVector, std::vector<PointPtr> *pointsVector)
{
    int dimension;
    for (int i = 0; i < linesVector.size(); i++)
    {
        // separate std::string by Spaces
        // pick every element from 2nd to std::endl
        // read point coordinates
        PointPtr currPoint = new Point;
        std::string word = "";
        dimension = 0;
        for (char x : linesVector[i])
        {
            if (x == ' ')
            {
                if (dimension)
                    currPoint->coords.push_back(atof(word.c_str()));
                else
                    currPoint->id = word;
                word = "";

                dimension++;
            }
            else
            {
                word = word + x;
            }
        }

        pointsVector->push_back(currPoint);
    }
    dimension--;
    return dimension;
}

void sort_neighbours(kNeighboursPtr k_nearest_neighbours, int k_neighbours) // sort distance in a vector of k distances
{
    int k = k_neighbours; // number of neighbours
    NeighbourPtr tempNeighbour;

    for (int i = k - 1; i > 0; i--)
    {
        if (k_nearest_neighbours->neighbours[i]->dist < k_nearest_neighbours->neighbours[i - 1]->dist)
        {
            tempNeighbour = k_nearest_neighbours->neighbours[i - 1];
            k_nearest_neighbours->neighbours[i - 1] = k_nearest_neighbours->neighbours[i];
            k_nearest_neighbours->neighbours[i] = tempNeighbour;
        }
    }
}

void sort_points(std::vector<PointPtr> *Data) // sort distance in a vector of k distances
{
    int k = Data->size();
    PointPtr tempPoint;
    for (int i = k - 1; i > 0; i--)
    {
        if ((*Data)[i]->id < (*Data)[i - 1]->id)
        {
            tempPoint = (*Data)[i];
            (*Data)[i] = (*Data)[i - 1];
            (*Data)[i - 1] = tempPoint;
        }
    }
}

void sort_points_str(std::vector<std::string> *Data)
{
    int k = Data->size();
    std::string tempPoint;
    for (int i = k - 1; i > 0; i--)
    {
        if ((*Data)[i] < (*Data)[i - 1])
        {
            tempPoint = (*Data)[i];
            (*Data)[i] = (*Data)[i - 1];
            (*Data)[i - 1] = tempPoint;
        }
    }
}

int notAlreadyExists(kNeighboursPtr k_nearest_neighbours, std::string pointID)
{

    for (int i = 0; i < k_nearest_neighbours->size; i++)
        if (k_nearest_neighbours->neighbours[i]->point->id == pointID)
            return 0;
    return 1;
}

kNeighboursPtr find_k_true_neighbours(PointPtr queryPoint, int k_neighbours, std::vector<PointPtr> inputPoints, int dim)
{
    NeighbourPtr currNeighbour = new Neighbour;
    kNeighboursPtr returnData = new kNeighbours;

    returnData->neighbours.resize(k_neighbours);
    returnData->size = k_neighbours;

    for (int i = 0; i < k_neighbours; i++)
    {
        returnData->neighbours[i] = new Neighbour;
        returnData->neighbours[i]->point = NULL;
        returnData->neighbours[i]->dist = INT32_MAX; // initialize distance with a very big value
    }

    for (int i = 0; i < inputPoints.size(); i++)
    {
        currNeighbour->point = inputPoints[i];
        currNeighbour->dist = euclideanDistance(queryPoint, currNeighbour->point, dim);

        if (currNeighbour->dist < returnData->neighbours[k_neighbours - 1]->dist && currNeighbour->dist > 0)
        {
            // if (returnData->size < k_neighbours)
            //     returnData->size++;
            returnData->neighbours[k_neighbours - 1]->point = currNeighbour->point;
            returnData->neighbours[k_neighbours - 1]->dist = currNeighbour->dist;

            if (k_neighbours > 1)
                sort_neighbours(returnData, k_neighbours);
        }
    }

    delete currNeighbour;
    return returnData;
}

std::string checkRerun()
{
    std::cout << "Rerun Program?..." << std::endl
              << "======Options======" << std::endl
              << "CONT to rerun" << std::endl
              << "TERM to terminate" << std::endl
              << "===================" << std::endl;
    std::string option;
    char ch;
    while (true)
    {
        option = "";
        while ((ch = getchar()) != '\n')
            option += ch;
        if (option == "CONT" || option == "TERM" || option == "cont" || option == "term")
            break;
    }
    return option;
}