#include "misc.h"

bool MetricCheck(std::string m)
{
	if (!m.compare("hamming") || !m.compare("lee")) return true;
	return false;
	//throw std::invalid_argument("This metric is not offered. Allowed metrics are hamming and lee.");
}

bool ParamCheck(double p)
{
    if (p >= 0 && p <= 1.0) return true;
    return false;
    //throw std::invalid_argument("Normalized value needs to be in the interval [0,1].");
}

bool AlgCheck(std::string a)
{
    if (!a.compare("prange") || !a.compare("dumer") || !a.compare("wagner")) return true;
    return false;
    //throw std::invalid_argument("This algorithm is not offered. Allowed algorithms are prange, dumer and wagner.");

}

bool AlphabetSizeCheck(int alphabetSize)
{
    if (alphabetSize >= 2) return true;
    return false;
    //throw std::invalid_argument("Alphabet size needs to be greater or equal to 2.");
}

double GoldenSectionSearch(double a, double b, double tol, const std::function<double(double)>& f)
{
    const double gr = (sqrt(5) + 1) / 2;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;
    while (abs(c - d) > tol)
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