#ifndef SPACE_H
#define SPACE_H

#include "fusion.h"
#include "monty.h"
#include "misc.h"

using namespace mosek::fusion;
using namespace monty;

class VectorSpace
{
protected:
	std::string metric;
	unsigned int alphabetSize;
public:
	VectorSpace(const std::string& m = "hamming", unsigned int as = 2)
	{
		if (AlphabetSizeCheck(as))
			alphabetSize = as;
		else
			throw std::invalid_argument("Alphabet size needs to be greater or equal to 2.");

		if(MetricCheck(m))
			metric = m;
		else
			throw std::invalid_argument("This metric is not offered. Allowed metrics are hamming and lee.");
	}
	VectorSpace(const VectorSpace& vs): metric(vs.metric), alphabetSize(vs.alphabetSize) {}
	VectorSpace& operator=(const VectorSpace& vs)
	{
		if (this == &vs)
			return *this;

		metric = vs.metric;
		alphabetSize = vs.alphabetSize;

		return *this;
	}
	std::string GetMetric() const { return metric; }
	void SetMetric(const std::string& m)
	{
		if(MetricCheck(m))
			metric = m;
		else
			throw std::invalid_argument("This metric is not offered. Allowed metrics are hamming and lee.");
	}
	unsigned int GetAlphabetSize() const { return alphabetSize; }
	void SetAlphabetSize(unsigned int as)
	{
		if(AlphabetSizeCheck(as))
			alphabetSize = as;
		else
			throw std::invalid_argument("Alphabet size needs to be greater or equal to 2.");
	}
	friend void LinConstrs(double, const VectorSpace&, const Model::t&, const Variable::t&);
	friend void ExpConicConstrs(const VectorSpace&, const Model::t&, const Variable::t&, const Variable::t&);
	double SphereSurfArea(double) const;
	friend double AvgVectorWeight(const VectorSpace&, double);
	friend unsigned int Weight(const VectorSpace&, unsigned int);
	friend unsigned int MaxWeight(const VectorSpace&);
};

#endif
