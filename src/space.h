#ifndef SPACE_H
#define SPACE_H

#include "monty.h"
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

void MetricCheck(std::string);

class VectorSpace
{
protected:
	std::string metric;
	unsigned int alphabetSize;
public:
	VectorSpace(const std::string& m = "hamming", unsigned int as = 2) : alphabetSize(as)
	{
		MetricCheck(m);
		metric = m;
	}
	VectorSpace(const VectorSpace& vs): metric(vs.metric)
	{
		alphabetSize = vs.alphabetSize;
	}
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
		MetricCheck(m);
		metric = m;
	}
	unsigned int GetAlphabetSize() const { return alphabetSize; }
	void SetAlphabetSize(unsigned int as) { alphabetSize = as; }
	friend void LinConstrs(double, const VectorSpace&, const Model::t&, const Variable::t&);
	friend void ExpConicConstrs(const VectorSpace&, const Model::t&, const Variable::t&, const Variable::t&);
	double SphereSurfArea(double) const;
	friend double AvgVectorWeight(const VectorSpace&, double);
	friend unsigned int Weight(const VectorSpace&, unsigned int);
	friend unsigned int MaxWeight(const VectorSpace&);
};

#endif
