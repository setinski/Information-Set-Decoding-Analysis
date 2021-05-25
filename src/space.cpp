#include <vector>
#include <algorithm>
#include "space.h"

void LinConstrs(double distance, const VectorSpace& space, const Model::t& M, const Variable::t& x)
{
	M->constraint("lc1", Expr::sum(x), Domain::equalsTo(1.0));

	std::vector<double> coeffs(space.alphabetSize);
	for (auto i = 0; i < space.alphabetSize; ++i)
	{
		coeffs[i] = (double) Weight(space, i);
	}

	auto coeffs_ptr = new_array_ptr<double>(coeffs);
	M->constraint("lc2", Expr::dot(coeffs_ptr, x), Domain::equalsTo(distance * MaxWeight(space)));
}

void ExpConicConstrs(const VectorSpace& space, const Model::t& M, const Variable::t& t, const Variable::t& x)
{
	for (auto i = 0; i < space.alphabetSize; i++)
	{
		M->constraint(Expr::hstack(1, x->index(i), t->index(i)), Domain::inPExpCone());
	}
}

double VectorSpace::SphereSurfArea(double distance) const
{
	Model::t M = new Model("space");
	auto _M = finally([&]() {M->dispose(); });

	Variable::t x = M->variable("x", alphabetSize, Domain::greaterThan(0.0));
	Variable::t t = M->variable("t", alphabetSize);

	LinConstrs(distance, *this, M, x);	
	ExpConicConstrs(*this, M, t, x);

	M->objective("obj", ObjectiveSense::Maximize, Expr::sum(t));
	M->solve();

	if (M->getProblemStatus() == ProblemStatus::PrimalAndDualFeasible)
	{
		if (M->getPrimalSolutionStatus() == SolutionStatus::Optimal)
		{
			return M->primalObjValue() / log(alphabetSize);
		}
	}
	return -1;
}

unsigned int Weight(const VectorSpace& space, unsigned int element)
{
	element %= space.alphabetSize;

	unsigned int weight;
	if (!space.metric.compare("hamming"))
	{
		weight = element ? 1 : 0;
	}
	else if (!space.metric.compare("lee"))
	{
		weight = std::min(element % space.alphabetSize, space.alphabetSize - element % space.alphabetSize);
	}
	else
	{
		throw std::invalid_argument("Invalid metric (needs to be either Hamming or Lee).");
	}
	return weight;
}

unsigned int MaxWeight(const VectorSpace& space)
{
	std::vector<unsigned int> weights(space.alphabetSize);
	for (auto i = 0; i < space.alphabetSize; ++i)
	{
		weights[i] = Weight(space, i);
	}

	return *std::max_element(weights.begin(), weights.end());
}

double AvgVectorWeight(const VectorSpace& space, double length)
{
	double avgWeight = 0;
	for (auto i = 0; i < space.alphabetSize; ++i)
	{
		avgWeight += (double)Weight(space, i);
	}

	return avgWeight / space.alphabetSize * length;
}