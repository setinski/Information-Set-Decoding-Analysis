#ifndef ISD_H
#define ISD_H

#include "space.h"
#include <string>
#include <functional>

extern const double tol;
extern const double epsilon;
extern bool quantum;

void ParamCheck(double);
void AlgCheck(std::string);
void AlphabetSizeCheck(int);

double AvgSolsNum(const VectorSpace&, double, double);
double UpperRoot(std::string, unsigned int, double);
double LowerRoot(std::string, unsigned int, double);
double RunTime(std::string, std::string, unsigned int, double,
				double&, double&, double&, unsigned int&);

double GoldenSectionSearch(double, double, double,
						const std::function<double(double)>&);

class InformationSetDecoding
{
private:
	int alphabetSize;
	double codeRate, weight, surfaceW;
	std::string metric, algorithm;
	VectorSpace space;
public:
	InformationSetDecoding(unsigned int as = 2, double cr = 0.5, double w = 0.5,
		const std::string& m = "hamming", const std::string& a = "prange");
	InformationSetDecoding(const InformationSetDecoding& isd) : 
		alphabetSize(isd.alphabetSize), codeRate(isd.codeRate), 
		weight(isd.weight), surfaceW(isd.surfaceW), metric(isd.metric),
		algorithm(isd.algorithm), space(isd.space){}
	InformationSetDecoding& operator=(const InformationSetDecoding&);
	unsigned int GetAlphaSize() const { return alphabetSize; }
	void SetAlphaSize(unsigned int as) { alphabetSize = as; }
	double GetCodeRate() const { return codeRate; }
	void SetCodeRate(double cr)
	{
		ParamCheck(cr);
		codeRate = cr;
	}
	double GetWeight() const { return weight; }
	void SetWeight(double w)
	{
		ParamCheck(w);
		weight = w;
	}
	std::string GetMetric() const { return metric; }
	void SetMetric(const std::string& m)
	{
		MetricCheck(m);
		metric = m;
	}
	std::string GetAlg() const { return algorithm; }
	void SetAlg(const std::string& a)
	{
		AlgCheck(a);
		algorithm = a;
	}
	double GetSurfaceW() const { return surfaceW; }
	VectorSpace GetSpace() const { return space; }
	double BDayDecCost(double, double, double) const;
	double SolsPerIter(double, double, double) const;
	double IterCost(double, double, double) const;
	double SolsNum() const;
	double PartSolProb(double, double) const;
	double AnySolProb(double, double) const;
	double RunTime(double, double, unsigned int&) const;
	double GoldenSectionSearch(double&, unsigned int&);
	double GoldenSectionSearch(double&, double&, unsigned int&);
	// add Destructor
};

#endif