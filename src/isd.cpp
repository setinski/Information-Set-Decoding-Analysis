#include "isd.h"
#include <algorithm>
#include <boost/math/tools/roots.hpp>

using boost::math::tools::bisect;

/* Constants: */

    /* Terminating tolerance. */
    const double tol = 1e-5;

    /* Tolerated computational error. */
    const double epsilon = 1e-5;

    /* Classical algorithm: quantum = false;
    Quantum algorithm: quantum = true. */
    bool quantum = false;
    //int bits = std::numeric_limits<double>::digits;

/* Parameters: in range [0,1]. */
void ParamCheck(double p)
{
    if (p >= 0 && p <= 1.0) return;

    throw std::invalid_argument("Normalized value needs to be in the interval [0,1].");

}

/* Algorithms compared: Prange's, Dumer/Stern's, Wagner's based. */
void AlgCheck(std::string a)
{
    if (!a.compare("prange") || !a.compare("dumer") || !a.compare("wagner")) return;

    throw std::invalid_argument("This algorithm is not offered. Allowed algorithms are prange, dumer and wagner.");

}

/* Alphabet size: greater or equal to 2. */
void AlphabetSizeCheck(int alphabetSize)
{
    if (alphabetSize >= 2) return;

    throw std::invalid_argument("Alphabet size needs to be greater or equal to 2.");
}

/* Instance of Information Set Decoding problem. */
InformationSetDecoding::InformationSetDecoding(unsigned int as, double cr, double w, const std::string& m, const std::string& a)
{
    alphabetSize = as;

    MetricCheck(m);
    metric = m;

    ParamCheck(cr);
    codeRate = cr;

    ParamCheck(w);
    weight = w;
    
    AlgCheck(a);
    algorithm = a;

    space = VectorSpace(metric, alphabetSize);
    surfaceW = space.SphereSurfArea(weight);
    if (surfaceW == -1)
    {
        throw std::invalid_argument("Surface area surfaceW is not found.");
    }
}

InformationSetDecoding& InformationSetDecoding::operator=(const InformationSetDecoding& isd)
{
    if (this == &isd)
        return *this;

    alphabetSize = isd.alphabetSize;
    codeRate = isd.codeRate;
    weight = isd.weight;
    metric = isd.metric;
    algorithm = isd.algorithm;
    space = isd.space;
    surfaceW = isd.surfaceW;

    return *this;
}

/* Cost of birthday decoding. */
double InformationSetDecoding::BDayDecCost(double paramL, double paramP, double optLevelNum) const
{
    double distance1 = paramP / (codeRate + paramL);
    try
    {
        ParamCheck(distance1);
    }
    catch (std::invalid_argument e)
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error("Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument("Surface area surface1 is not found.");
    }

    if (quantum)
    {
        return std::min(surface1 / (pow(2, optLevelNum) + 1), paramL / optLevelNum);
    }
    else
    {
        return std::min(surface1 / pow(2, optLevelNum), paramL / optLevelNum);
    }
}

/* Number of solutions per iteration. */
double InformationSetDecoding::SolsPerIter(double paramL, double paramP, double optLevelNum) const
{
    double distance1 = paramP / (codeRate + paramL);
    try
    {
        ParamCheck(distance1);
    }
    catch (std::invalid_argument e)
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error("Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument("Surface area surface1 is not found.");
    }

    double bottomListSize;
    if (quantum)
    {
        bottomListSize = std::min(surface1 / (pow(2, optLevelNum) + 1), paramL / optLevelNum);
    }
    else
    {
        bottomListSize = std::min(surface1 / pow(2, optLevelNum), paramL / optLevelNum);
    }

    double paramM = paramL - (optLevelNum - 1) * bottomListSize;
    if (quantum)
    {
        return 3 * bottomListSize - paramM;
    }
    else
    {
        return 2 * bottomListSize - paramM;
    }   
}

/* Iteration cost. */
double InformationSetDecoding::IterCost(double paramL, double paramP, double optLevelNum) const
{   
    double distance1 = paramP / (codeRate + paramL);
    try
    {
        ParamCheck(distance1);
    }
    catch (std::invalid_argument e)
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error("Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument("Surface area surface1 is not found.");
    } 

    if (quantum)
    {
        return std::min(surface1 / (pow(2, optLevelNum) + 1), paramL / optLevelNum);
    }
    else
    {
        return std::min(surface1 / pow(2, optLevelNum), paramL / optLevelNum);
    }
}

/* Expected number of solutions. */
double InformationSetDecoding::SolsNum() const
{
    return std::max(surfaceW - (1 - codeRate), 0.0);
}

/* Probability of finding a particular solution. */
double InformationSetDecoding::PartSolProb(double paramL, double paramP) const
{
    double distance2 = (weight - paramP) / (1 - codeRate - paramL);
    try
    {
        ParamCheck(distance2);
    }
    catch (std::invalid_argument e)
    {
        if (distance2 >= 1.0 && distance2 <= 1.0 + epsilon)
        {
            distance2 = 1.0;
        }
        else
        {
            std::cout << distance2 << std::endl;
            throw std::runtime_error("Numerical issue on parameter (weight - paramP) / (1 - codeRate - paramL).");
        }
    }

    double surface2 = (1 - codeRate - paramL) * space.SphereSurfArea(distance2);
    if (surface2 == -1)
    {
        throw std::invalid_argument("Surface area surface2 is not found.");
    }

    return surface2 - surfaceW + paramL;
}

/* Probablility of findiing any solution. */
double InformationSetDecoding::AnySolProb(double paramL, double paramP) const
{
    double distance2 = (weight - paramP) / (1 - codeRate - paramL);
    try
    {
        ParamCheck(distance2);
    }
    catch (std::invalid_argument e)
    {
        if (distance2 >= 1.0 && distance2 <= 1.0 + epsilon)
        {
            distance2 = 1.0;
        }
        else
        {
            std::cout << distance2 << std::endl;
            throw std::runtime_error("Numerical issue on parameter (weight - paramP) / (1 - codeRate - paramL).");
        }
    }

    double surface2 = (1 - codeRate - paramL) * space.SphereSurfArea(distance2);
    if (surface2 == -1)
    {
        throw std::invalid_argument("Surface area surface2 is not found.");
    }

    return std::min(0.0, surface2 + paramL - std::min(1 - codeRate, surfaceW));
}

/* Running time of the algorithm. */
double InformationSetDecoding::RunTime(double paramL, double paramP, unsigned int& optLevelNum) const
{
    double distance1 = paramP / (codeRate + paramL);
    try
    {
        ParamCheck(distance1);
    }
    catch (std::invalid_argument e)
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error("Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument("Surface area surface1 is not found.");
    }

    double distance2 = (weight - paramP) / (1 - codeRate - paramL);
    try
    {
        ParamCheck(distance2);
    }
    catch (std::invalid_argument e)
    {
        if (distance2 >= 1.0 && distance2 <= 1.0 + epsilon)
        {
            distance2 = 1.0;
        }
        else
        {
            std::cout << distance2 << std::endl;
            throw std::runtime_error("Numerical issue on parameter (weight - paramP) / (1 - codeRate - paramL).");
        }
    }

    double surface2 = (1 - codeRate - paramL) * space.SphereSurfArea(distance2);
    if (surface2 == -1)
    {
        throw std::invalid_argument("Surface area surface2 is not found.");
    }

    auto theorem = [=](unsigned int optLevelNum) -> bool
    {
        return bool(paramL <= optLevelNum / pow(2.0, optLevelNum) * surface1);
    };

    if (!algorithm.compare("prange") || !algorithm.compare("dumer"))
    {
        optLevelNum = 1;
    }
    else
    {
        optLevelNum = 2;
        while (theorem(optLevelNum))
        {
            optLevelNum++;
        }
        optLevelNum--;
    }

    double bottomListSize;
    if (quantum)
    {
        bottomListSize = std::min(surface1 / (pow(2, optLevelNum) + 1), paramL / optLevelNum);
    }
    else
    {
        bottomListSize = std::min(surface1 / pow(2, optLevelNum), paramL / optLevelNum);
    }

    double paramM = paramL - (optLevelNum - 1) * bottomListSize;
    if (quantum)
    {
        return bottomListSize - 0.5 * std::min(0.0, std::min(0.0, surface2 + paramL - std::min(1 - codeRate, surfaceW)) + 3 * bottomListSize - paramM);
    }
    else
    {
        return bottomListSize - std::min(0.0, std::min(0.0, surface2 + paramL - std::min(1 - codeRate, surfaceW)) + 2 * bottomListSize - paramM);
    }
}

/* Golden section search. */
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

/* Golden section search on RunTime function with one parameter
(used for Prange's algorithm parameter optimization). */
double InformationSetDecoding::GoldenSectionSearch(double& paramP, unsigned int& optLevelNum)
{
    double paramPLow = std::max(0.0, weight - (1 - codeRate));
    double paramPHigh = std::min(weight, codeRate);
    std::function<double(double)> runTimeFun = [&](double p) {return RunTime(0.0, p, optLevelNum); };
    paramP = ::GoldenSectionSearch(paramPLow, paramPHigh, tol, runTimeFun);
    
    return log2(alphabetSize) * RunTime(0.0, paramP, optLevelNum);
}

/* Golden section search on RunTime function with two parameters
(used for Dumer/Stern's and Wagner's algorithm. )*/
double InformationSetDecoding::GoldenSectionSearch(double& paramL, double& paramP, unsigned int& optLevelNum)
{
    //const double paramLLow = 0, paramLHigh = std::min(1 - codeRate - 0.1, 0.25);
    const double paramLLow = 0, paramLHigh = 1 - codeRate;
    double a = paramLLow, b = paramLHigh;
    const double gr = (sqrt(5) + 1) / 2;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;

    while (abs(c - d) > tol)
    {
        double paramPLowC = std::max(0.0, weight - (1 - codeRate - c));
        double paramPHighC = std::min(weight, codeRate + c);
        std::function<double(double)> runTimeC = [&](double p) {return RunTime(c, p, optLevelNum); };
        double optC = ::GoldenSectionSearch(paramPLowC, paramPHighC, tol, runTimeC);

        double paramPLowD = std::max(0.0, weight - (1 - codeRate - d));
        double paramPHighD = std::min(weight, codeRate + d);
        std::function<double(double)> runTimeD = [&](double p) {return RunTime(d, p, optLevelNum); };
        double optD = ::GoldenSectionSearch(paramPLowD, paramPHighD, tol, runTimeD);

        if (RunTime(c, optC, optLevelNum) < RunTime(d, optD, optLevelNum))
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

    paramL = (b + a) / 2;
    std::function<double(double)> runTime = [&](double p) {return RunTime(paramL, p, optLevelNum); };
    paramP = ::GoldenSectionSearch(std::max(0.0, weight - (1 - codeRate - paramL)),
        std::min(weight, codeRate + paramL), tol, runTime);

    return log2(alphabetSize) * RunTime(paramL, paramP, optLevelNum);
}

/* Termination condition. */
struct TerminationCondition {
    bool operator() (double min, double max) {
        return abs(min - max) <= tol;
    }
};

/* Average number of solutions per iteration of the algorithm. */
double AvgSolsNum(const VectorSpace& space, double weight, double codeRate)
{
    return space.SphereSurfArea(weight) - (1 - codeRate);
}

/* Calculating the upper root of the average number of solutions given as a function on weight. */
double UpperRoot(std::string metric, unsigned int alphabetSize, double codeRate)
{
    VectorSpace space(metric, alphabetSize);
    std::function<double(double)> avgSolsNumCR = [&](double cr) {return AvgSolsNum(space, 1 - epsilon, cr); };
    std::pair<double, double> bracketsCR = bisect(avgSolsNumCR, 0.0, 1.0, TerminationCondition());
    double highestCodeRate = (bracketsCR.first + bracketsCR.second) / 2;

    if (codeRate > highestCodeRate)
        return 1;
    else
    {
        std::function<double(double)> avgSolsNumW = [&](double w) {return AvgSolsNum(space, w, codeRate); };
        std::pair<double, double> bracketsW = bisect(avgSolsNumW, AvgVectorWeight(space, 1.0) / MaxWeight(space), 1.0, TerminationCondition());
        return (bracketsW.first + bracketsW.second) / 2;
    }
}

/* Calculating the lower root of the average number of solutions given as a function on weight. */
double LowerRoot(std::string metric, unsigned int alphabetSize, double codeRate)
{
    VectorSpace space(metric, alphabetSize);
    std::function<double(double)> avgSolsNumW = [&](double w) {return AvgSolsNum(space, w, codeRate); };
    std::pair<double, double> bracketsW = bisect(avgSolsNumW, 0.0, AvgVectorWeight(space, 1.0) / MaxWeight(space), TerminationCondition());
    return (bracketsW.first + bracketsW.second) / 2;
}

/* Running time of the algorithm that also gives the associated weight and corresponding paramL, paramP and optLevelNum. */
double RunTime(std::string metric, std::string algorithm, unsigned int alphabetSize, double codeRate,
    double& weight, double& paramL, double& paramP, unsigned int& optLevelNum)
{
    double runTime;
    if (!metric.compare("hamming") && alphabetSize != 3)
    {
        weight = LowerRoot(metric, alphabetSize, codeRate) - epsilon;
    }
    else
    {
        weight = UpperRoot(metric, alphabetSize, codeRate) - epsilon;
    }

    InformationSetDecoding isd(alphabetSize, codeRate, weight, metric, algorithm);
    if (!algorithm.compare("prange"))
    {
        try
        {
            paramL = 0;
            runTime = isd.GoldenSectionSearch(paramP, optLevelNum);
        }
        catch (std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
            return -1;
        }
    }
    else if (!algorithm.compare("dumer") || !algorithm.compare("wagner"))
    {
        try
        {
            runTime = isd.GoldenSectionSearch(paramL, paramP, optLevelNum);
        }
        catch (std::runtime_error& e)
        {
            std::cout << e.what() << std::endl;
            return -1;
        }
    }
    else
    {
        std::cout << "This algorithm is not offered. Allowed algorithms are prange, dumer and wagner." << std::endl;
        return -1;
    }
    return -runTime;
}
