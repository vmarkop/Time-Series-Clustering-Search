#include "mathUtils.h"

std::default_random_engine generator;

double normalDistributionGenerator(const double mi, const double sigma)
{
    // Source: https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
    std::normal_distribution<double> distribution(mi, sigma);
    // cout << "NormalDistGen " << distribution(generator) << endl;
    return distribution(generator);
}

double euclideanDistance(const PointPtr x, const PointPtr y, int dimension)
{
    double dist = 0.0;
    for (int i = 0; i < dimension; i++)
        dist += (x->coords[i] - y->coords[i]) * (x->coords[i] - y->coords[i]);
    return sqrt(dist);
}

double DFDistance(PointPtr x, PointPtr y, int dimension)
{
}

double uniformDistributionGenerator(const double alpha, const double beta)
{
    // Source: https://www.cplusplus.com/reference/random/normal_distribution/
    std::random_device randomDevice;
    std::mt19937 generator(randomDevice());
    std::uniform_real_distribution<> distribution(alpha, beta);
    //cout << "UniformDistGen " << distribution(generator) << endl;

    return distribution(generator);
}

int avoidOverFlowModulo(int a, int b, int m, char op)
{
    switch (op)
    {
    case '+':
        return euclideanModulo(euclideanModulo(a, m) + euclideanModulo(b, m), m);
    case '-':
        return euclideanModulo(euclideanModulo(a, m) - euclideanModulo(b, m), m);
    case '*':
        return euclideanModulo(euclideanModulo(a, m) * euclideanModulo(b, m), m);
    case '/':
        if (b == 0)
        {
            std::cerr << "Division by zero" << std::endl;
            exit(1);
        }
        return euclideanModulo(euclideanModulo(a, m) / euclideanModulo(b, m), m);
    case '%':
        return euclideanModulo(euclideanModulo(a, m) % euclideanModulo(b, m), m);
        ;
    default:
        std::cerr << "Wrong Operator" << std::endl;
        exit(1);
    }
}

int euclideanModulo(int x, int y)
{
    int returnValue = x % y;
    return returnValue >= 0 ? returnValue : returnValue + std::abs(y);
}

int powerWithBase2(int exp)
{

    int retValue = 1;

    for (int i = 0; i < exp; i++)
    {
        retValue *= 2;
    }

    return retValue;
}

bool is_number(const std::string &str)
// https://stackoverflow.com/questions/4654636/how-to-determine-if-a-string-is-a-number-with-c
{
    std::string::const_iterator iteration = str.begin();
    while (iteration != str.end() && std::isdigit(*iteration))
        ++iteration;
    return !str.empty() && iteration == str.end();
}