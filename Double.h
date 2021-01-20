//
//  asinh, acosh, atanh - C++11 fuctions that are not in math.h
//  D.P. Mitchell  2018/09/20.
//
#include "stdafx.h"

double
asinh(double x)
{
    return log(x + sqrt(x*x + 1.0));
}

double
acosh(double x)
{
    return log(x + sqrt(x*x - 1.0));
}

double
atanh(double x)
{
    return 0.5*log((1.0 + x) / (1.0 - x));
}