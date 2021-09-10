#include "misc.h"
#include <cmath>

bool MetricCheck(std::string m)
{
	if (!m.compare("hamming") || !m.compare("lee")) return true;
	return false;
}

bool ParamCheck(double p)
{
    if (p >= 0 && p <= 1.0) return true;
    return false;
}

bool AlgCheck(std::string a)
{
    if (!a.compare("prange") || !a.compare("dumer") || !a.compare("wagner")) return true;
    return false;
}

bool AlphabetSizeCheck(int alphabetSize)
{
    if (alphabetSize >= 2) return true;
    return false;
}

double GoldenSectionSearch(double a, double b, double tol, const std::function<double(double)>& f)
{
    const double gr = (sqrt(5) + 1) / 2;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;
    while (std::abs(c - d) > tol)
    {
        if (f(c) < f(d))
        {
            b = d;
        }
        else
        {
            a = c;
        }
        c = b - (b - a) / gr;
        d = a + (b - a) / gr;
    }

    return (b + a) / 2;
}
