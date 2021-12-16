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

double min(double x1, double x2, double x3)
{
    double retMin = x1;
    if (x2 < retMin)
    {
        retMin = x2;
    }
    if (x3 < retMin)
    {
        retMin = x3;
    }
    return retMin;
}

int minIndx(double x1, double x2, double x3)
{
    double retMin = x1;
    int indxMin = 0;
    if (x2 < retMin)
    {
        retMin = x2;
        indxMin = 1;
    }
    if (x3 < retMin)
    {
        retMin = x3;
        indxMin = 2;
    }
    return indxMin;
}

double DFDistance(PointPtr p, PointPtr q, int dimension, std::vector<std::vector<double>> *_c)
{
    if (_c == NULL)
        _c = new std::vector<std::vector<double>>;
    std::vector<PointPtr> p_points, q_points;
    p_points.resize(dimension / 2);
    q_points.resize(dimension / 2);
    for (int i = 0; i < dimension / 2; i++)
    {
        p_points[i] = new PointStruct;
        p_points[i]->id = p->id;
        p_points[i]->coords.resize(2);
        p_points[i]->coords[0] = p->coords[i * 2];
        p_points[i]->coords[1] = p->coords[(i * 2) + 1];

        q_points[i] = new PointStruct;
        q_points[i]->id = q->id;
        q_points[i]->coords.resize(2);
        q_points[i]->coords[0] = q->coords[i * 2];
        q_points[i]->coords[1] = q->coords[(i * 2) + 1];
    }

    // Initializing _c dimention
    _c->resize(dimension / 2);
    for (int i = 0; i < dimension / 2; i++)
    {
        (*_c)[i].resize(dimension / 2);
    }
    (*_c)[0][0] = euclideanDistance(p_points[0], q_points[0], 2);
    for (int i = 0; i < dimension / 2; i++)
    {
        for (int j = i; j < dimension / 2; j++)
        {
            if (i == 0 && j > 0)
            {
                // Creating row
                if ((*_c)[0][j - 1] > euclideanDistance(p_points[0], q_points[j], 2))
                    (*_c)[0][j] = (*_c)[0][j - 1];
                else
                    (*_c)[0][j] = euclideanDistance(p_points[0], q_points[j], 2);

                // Creating column
                if (i != j)
                {
                    int temp = i;
                    i = j;
                    j = temp;

                    if ((*_c)[i - 1][0] > euclideanDistance(p_points[i], q_points[0], 2))
                        (*_c)[i][0] = (*_c)[i - 1][0];
                    else
                        (*_c)[i][0] = euclideanDistance(p_points[i], q_points[0], 2);
                    temp = i;
                    i = j;
                    j = temp;
                }
            }
            else if (i > 0 && j > 0)
            {
                // Creating row
                double min_c = min((*_c)[i][j - 1], (*_c)[i - 1][j], (*_c)[i - 1][j - 1]);
                if (min_c > euclideanDistance(p_points[i], q_points[j], 2))
                    (*_c)[i][j] = min_c;
                else
                    (*_c)[i][j] = euclideanDistance(p_points[i], q_points[j], 2);

                // Creating column
                if (i != j)
                {
                    int temp = i;
                    i = j;
                    j = temp;
                    min_c = min((*_c)[i][j - 1], (*_c)[i - 1][j], (*_c)[i - 1][j - 1]);
                    if (min_c > euclideanDistance(p_points[i], q_points[j], 2))
                        (*_c)[i][j] = min_c;
                    else
                        (*_c)[i][j] = euclideanDistance(p_points[i], q_points[j], 2);
                    temp = i;
                    i = j;
                    j = temp;
                }
            }
        }
    }
    return (*_c)[(dimension / 2) - 1][(dimension / 2) - 1];
}

double uniformDistributionGenerator(const double alpha, const double beta)
{
    // Source: https://www.cplusplus.com/reference/random/normal_distribution/
    std::random_device randomDevice;
    std::mt19937 generator(randomDevice());
    std::uniform_real_distribution<> distribution(alpha, beta);
    // cout << "UniformDistGen " << distribution(generator) << endl;

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