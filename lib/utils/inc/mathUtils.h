#ifndef _MATH_UTILS_H_
#define _MATH_UTILS_H_

#include <cmath>
#include <iostream>
#include <random>
#include <string.h>
#include "projectUtils.h"

double normalDistributionGenerator(double mi = 0.0, double sigma = 1.0);
double uniformDistributionGenerator(double alpha = 0.0, double beta = 1.0);
double euclideanDistance(PointPtr x, PointPtr y, int dimension);
double min(double x1, double x2, double x3);
int minIndx(double x1, double x2, double x3);
double DFDistance(std::vector<std::vector<double>> *_c, PointPtr p, PointPtr q, int dimension);
int avoidOverFlowModulo(int a, int b, int m, char op);
int euclideanModulo(int x, int y);
int powerWithBase2(int exp);
bool is_number(const std::string &s);

#endif