#include "isd.h"
#include <algorithm>
#include <boost/math/tools/roots.hpp>

using boost::math::tools::bisect;

/*---------------------------------ISD functions------------------------------------*/

InformationSetDecoding::InformationSetDecoding(unsigned int as, double cr, double w,
                                         const std::string& m, const std::string& a)
{
    if(AlphabetSizeCheck(as))
        alphabetSize = as;
    else
        throw std::invalid_argument(
            "Invalid aphabet size (needs to be an integer value greater than 2).");

    if(MetricCheck(m))
        metric = m;
    else
        throw std::invalid_argument(
            "Invalid metric (needs to be either Hamming or Lee).");

    if(ParamCheck(cr))
        codeRate = cr;
    else
        throw std::invalid_argument(
            "Normalized value needs to be in the interval [0,1].");
    
    if(ParamCheck(w))
        weight = w;
    else
        throw std::invalid_argument(
            "Normalized value needs to be in the interval [0,1].");
    
    if(AlgCheck(a))
        algorithm = a;
    else
        throw std::invalid_argument(
            "This algorithm is not offered. Allowed algorithms are prange, dumer and wagner.");

    space = VectorSpace(metric, alphabetSize);
    surfaceW = space.SphereSurfArea(weight);
    if (surfaceW == -1)
    {
        throw std::invalid_argument(
            "Surface area surfaceW is not found.");
    }
}

InformationSetDecoding& InformationSetDecoding::operator=
                        (const InformationSetDecoding& isd)
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

double InformationSetDecoding::BDayDecCost(double paramL, double paramP,
                                            double optLevelNum) const
{
    double distance1 = paramP / (codeRate + paramL);
    if(!ParamCheck(distance1))
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error(
                "Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument(
            "Surface area surface1 is not found.");
    }

    if (quantum)
    {
        return std::min(surface1 / (pow(2, optLevelNum) + 1),
                        paramL / optLevelNum);
    }
    else
    {
        return std::min(surface1 / pow(2, optLevelNum),
                        paramL / optLevelNum);
    }
}

double InformationSetDecoding::SolsPerIter(double paramL, double paramP,
                                            double optLevelNum) const
{
    double distance1 = paramP / (codeRate + paramL);
    
    if(!ParamCheck(distance1))
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error(
                "Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument(
            "Surface area surface1 is not found.");
    }

    double bottomListSize;
    if (quantum)
    {
        bottomListSize = std::min(surface1 / (pow(2, optLevelNum) + 1),
                                    paramL / optLevelNum);
    }
    else
    {
        bottomListSize = std::min(surface1 / pow(2, optLevelNum),
                                    paramL / optLevelNum);
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

double InformationSetDecoding::IterCost(double paramL, double paramP,
                                        double optLevelNum) const
{   
    double distance1 = paramP / (codeRate + paramL);
    if(!ParamCheck(distance1))
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error(
                "Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument(
            "Surface area surface1 is not found.");
    } 

    if (quantum)
    {
        return std::min(surface1 / (pow(2, optLevelNum) + 1),
                        paramL / optLevelNum);
    }
    else
    {
        return std::min(surface1 / pow(2, optLevelNum),
                        paramL / optLevelNum);
    }
}

double InformationSetDecoding::SolsNum() const
{
    return std::max(surfaceW - (1 - codeRate), 0.0);
}

double InformationSetDecoding::PartSolProb(double paramL, double paramP) const
{
    double distance2 = (weight - paramP) / (1 - codeRate - paramL);
    if(!ParamCheck(distance2))
    {
        if (distance2 >= 1.0 && distance2 <= 1.0 + epsilon)
        {
            distance2 = 1.0;
        }
        else
        {
            std::cout << distance2 << std::endl;
            throw std::runtime_error(
                "Numerical issue on parameter (weight - paramP) / (1 - codeRate - paramL).");
        }
    }

    double surface2 = (1 - codeRate - paramL) * space.SphereSurfArea(distance2);
    if (surface2 == -1)
    {
        throw std::invalid_argument("Surface area surface2 is not found.");
    }

    return surface2 - surfaceW + paramL;
}

double InformationSetDecoding::AnySolProb(double paramL, double paramP) const
{
    double distance2 = (weight - paramP) / (1 - codeRate - paramL);
    if(!ParamCheck(distance2))
    {
        if (distance2 >= 1.0 && distance2 <= 1.0 + epsilon)
        {
            distance2 = 1.0;
        }
        else
        {
            std::cout << distance2 << std::endl;
            throw std::runtime_error(
                "Numerical issue on parameter (weight - paramP) / (1 - codeRate - paramL).");
        }
    }

    double surface2 = (1 - codeRate - paramL) * space.SphereSurfArea(distance2);
    if (surface2 == -1)
    {
        throw std::invalid_argument(
            "Surface area surface2 is not found.");
    }

    return std::min(0.0, surface2 + paramL - std::min(1 - codeRate, surfaceW));
}

double InformationSetDecoding::RunTime(double paramL, double paramP,
                                       unsigned int& optLevelNum) const
{
    double distance1 = paramP / (codeRate + paramL);
    if(!ParamCheck(distance1))
    {
        if (distance1 >= 1.0 && distance1 <= 1.0 + epsilon)
        {
            distance1 = 1.0;
        }
        else
        {
            std::cout << distance1 << std::endl;
            throw std::runtime_error(
                "Numerical issue on parameter paramP / (codeRate + paramL).");
        }
    }

    double surface1 = (codeRate + paramL) * space.SphereSurfArea(distance1);
    if (surface1 == -1)
    {
        throw std::invalid_argument(
            "Surface area surface1 is not found.");
    }

    double distance2 = (weight - paramP) / (1 - codeRate - paramL);
    if(!ParamCheck(distance2))
    {
        if (distance2 >= 1.0 && distance2 <= 1.0 + epsilon)
        {
            distance2 = 1.0;
        }
        else
        {
            std::cout << distance2 << std::endl;
            throw std::runtime_error(
                "Numerical issue on parameter (weight - paramP) / (1 - codeRate - paramL).");
        }
    }

    double surface2 = (1 - codeRate - paramL) * space.SphereSurfArea(distance2);
    if (surface2 == -1)
    {
        throw std::invalid_argument(
            "Surface area surface2 is not found.");
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
        bottomListSize = std::min(surface1 / (pow(2, optLevelNum) + 1),
                                  paramL / optLevelNum);
    }
    else
    {
        bottomListSize = std::min(surface1 / pow(2, optLevelNum),
                                  paramL / optLevelNum);
    }

    double paramM = paramL - (optLevelNum - 1) * bottomListSize;
    if (quantum)
    {
        return bottomListSize - 0.5 * std::min(0.0, std::min(0.0, surface2 + paramL
            - std::min(1 - codeRate, surfaceW)) + 3 * bottomListSize - paramM);
    }
    else
    {
        return bottomListSize - std::min(0.0, std::min(0.0, surface2 + paramL
            - std::min(1 - codeRate, surfaceW)) + 2 * bottomListSize - paramM);
    }
}

double InformationSetDecoding::GoldenSectionSearch(double& paramP,
                                         unsigned int& optLevelNum)
{
    double paramPLow = std::max(0.0, weight - (1 - codeRate));
    double paramPHigh = std::min(weight, codeRate);
    std::function<double(double)> runTimeFun = [&](double p)
                            {return RunTime(0.0, p, optLevelNum); };
    paramP = ::GoldenSectionSearch(paramPLow, paramPHigh, tol, runTimeFun);
    
    return log2(alphabetSize) * RunTime(0.0, paramP, optLevelNum);
}

double InformationSetDecoding::GoldenSectionSearch(double& paramL, double& paramP,
                                                   unsigned int& optLevelNum)
{
    const double gr = (sqrt(5) + 1) / 2;
    const double paramLLow = 0, paramLHigh = 1 - codeRate;
    double a = paramLLow, b = paramLHigh;
    double c = b - (b - a) / gr;
    double d = a + (b - a) / gr;

    while (abs(c - d) > tol)
    {
        double paramPLowC = std::max(0.0, weight - (1 - codeRate - c));
        double paramPHighC = std::min(weight, codeRate + c);
        std::function<double(double)> runTimeC = [&](double p)
                                    {return RunTime(c, p, optLevelNum); };
        double optC = ::GoldenSectionSearch(paramPLowC, paramPHighC, tol, runTimeC);

        double paramPLowD = std::max(0.0, weight - (1 - codeRate - d));
        double paramPHighD = std::min(weight, codeRate + d);
        std::function<double(double)> runTimeD = [&](double p)
                                    {return RunTime(d, p, optLevelNum); };
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
    std::function<double(double)> runTime = [&](double p)
                                    {return RunTime(paramL, p, optLevelNum); };
    paramP = ::GoldenSectionSearch(std::max(0.0, weight - (1 - codeRate - paramL)),
        std::min(weight, codeRate + paramL), tol, runTime);

    return log2(alphabetSize) * RunTime(paramL, paramP, optLevelNum);
}

/*-----------------------Hardest instances functions--------------------------------*/

struct TerminationCondition {
    bool operator() (double min, double max) {
        return abs(min - max) <= tol;
    }
};

double AvgSolsNum(const VectorSpace& space, double distance, double codeRate)
{
    double surface = space.SphereSurfArea(distance);
    if (surface == -1)
    {
        throw std::invalid_argument(
            "Surface area surface1 is not found.");
    }
    return space.SphereSurfArea(distance) - (1 - codeRate);
}

double UpperRoot(std::string metric, unsigned int alphabetSize, double codeRate)
{
    VectorSpace space(metric, alphabetSize);
    std::function<double(double)> avgSolsNumCR = [&](double cr)
                                    {return AvgSolsNum(space, 1 - epsilon, cr); };
    std::pair<double, double> bracketsCR = bisect(avgSolsNumCR, 0.0, 1.0,
                                    TerminationCondition());
    double highestCodeRate = (bracketsCR.first + bracketsCR.second) / 2;

    if (codeRate > highestCodeRate)
        return 1;
    else
    {
        std::function<double(double)> avgSolsNumW = [&](double w)
                                     {return AvgSolsNum(space, w, codeRate); };
        std::pair<double, double> bracketsW = bisect(avgSolsNumW,
            AvgVectorWeight(space, 1.0) / MaxWeight(space), 1.0, TerminationCondition());
        return (bracketsW.first + bracketsW.second) / 2;
    }
}

double LowerRoot(std::string metric, unsigned int alphabetSize, double codeRate)
{
    VectorSpace space(metric, alphabetSize);
    std::function<double(double)> avgSolsNumW = [&](double w)
                                        {return AvgSolsNum(space, w, codeRate); };
    std::pair<double, double> bracketsW = bisect(avgSolsNumW,
            0.0, AvgVectorWeight(space, 1.0) / MaxWeight(space), TerminationCondition());
    return (bracketsW.first + bracketsW.second) / 2;
}

double RunTime(std::string metric, std::string algorithm, unsigned int alphabetSize,
               double codeRate, double& weight, double& paramL, double& paramP,
               unsigned int& optLevelNum)
{
    double runTime;
    InformationSetDecoding isd(alphabetSize, codeRate, weight, metric, algorithm);
    
    if (!metric.compare("hamming") && alphabetSize != 3)
    {
        weight = LowerRoot(metric, alphabetSize, codeRate) - epsilon;
    }
    else
    {
        weight = UpperRoot(metric, alphabetSize, codeRate) - epsilon;
    }

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
        std::cout << "This algorithm is not offered. ";
        std::cout << "Allowed algorithms are prange, dumerand wagner.";
        std::cout << std::endl;
        return -1;
    }
    return -runTime;
}

/*--------------------------------------------------------------------------------*/
