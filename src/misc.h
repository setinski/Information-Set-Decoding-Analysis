#ifndef MISC_H
#define MISC_H
#include <string>
#include <functional>

bool MetricCheck(std::string);
bool ParamCheck(double);
bool AlgCheck(std::string);
bool AlphabetSizeCheck(int);

double GoldenSectionSearch(double, double, double, const std::function<double(double)>&);

#endif
