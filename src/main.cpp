#pragma omp parallel for
#pragma optimize( "", off )

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include "space.h"
#include "isd.h"

/* The command line parameters are: metric as a string ("lee" or "hamming"), algorithm as a string 
("prange", "dumer", "wagner") and if the algorithms is quantum as bool (0: classical, 1: quantum). */
void ArgcCheck(int argc)
{
	if (argc >= 3) return;

	throw std::invalid_argument("The number of arguments cannot be smaller than 3.");
}

void NumpointsCheck(int np)
{
	if (np >= 1) return;

	throw std::invalid_argument("Number of points needs to be greater than 1.");
}

int main(int argc, char** argv)
{
	// 1. Parameters' checks:
	try
	{
		ArgcCheck(argc);
	}
	catch (std::invalid_argument& e)
	{
		std::cerr << e.what() << std::endl;
		return -1;
	}

	std::string metric(argv[1]);
	try
	{
		MetricCheck(metric);
	}
	catch (std::invalid_argument& e)
	{
		std::cerr << e.what() << std::endl;
		return -1;
	}
	std::cout << "Metric:           " << metric << std::endl;

	std::string algorithm(argv[2]);
	try
	{
		AlgCheck(algorithm);
	}
	catch (std::invalid_argument& e)
	{
		std::cerr << e.what() << std::endl;
		return -1;
	}
	std::cout << "Algorithm:        " << algorithm << std::endl << std::endl;

	if (argc == 4)
	{
		quantum = atoi(argv[3]);
	}

	// 2. Calculating the hardest instances
	std::vector<unsigned int> alphabetSizes = { 3, 5, 7, 13, 23, 31, 43, 83, 163, 331, 643};
	std::string outputFilename = "results";
	if (quantum)
	{
		outputFilename += "_quantum";
	}
	std::fstream outputFile(outputFilename + ".txt", std::ofstream::out | std::ofstream::trunc);
	outputFile << "alphabetSize codeRate weight optLevelNum paramL paramP runtime(log 2) runtime(log alphabetSize)" << std::endl;
	outputFile.close();

	for (auto i = 0; i < alphabetSizes.size(); i++)
	{
		std::cout << "Alphabet size: " << alphabetSizes[i] << std::endl;
		auto t1 = std::chrono::high_resolution_clock::now();

		double paramL, paramP, weight;
		unsigned int optLevelNum;
		std::function<double(double)> runTimeFun = [&](double cr) {return RunTime(metric, algorithm, alphabetSizes[i],
																			cr, weight, paramL, paramP, optLevelNum); };
		double codeRate = GoldenSectionSearch(epsilon, 1.0 - epsilon, tol, runTimeFun);
		double runTime = RunTime(metric, algorithm, alphabetSizes[i], codeRate, weight, paramL, paramP, optLevelNum);

		outputFile.open(outputFilename + ".txt", std::ios_base::app);
		outputFile << std::fixed << std::setprecision(1) << std::setfill('0') << std::setw(3) << alphabetSizes[i] << "          ";
		outputFile << std::fixed << std::setprecision(3) << codeRate << "    ";
		outputFile << std::fixed << std::setprecision(3) << weight << "  ";
		outputFile << std::fixed << std::setprecision(1) << std::setfill('0') << std::setw(3) << optLevelNum << "         ";
		outputFile << std::fixed << std::setprecision(3) << paramL << "  ";
		outputFile << std::fixed << std::setprecision(3) << paramP << "  ";
		outputFile << std::fixed << std::setprecision(3) << -runTime << "                     ";
		outputFile << std::fixed << std::setprecision(3) << -runTime / log2(alphabetSizes[i]) << std::endl;
		outputFile.close();

		auto t2 = std::chrono::high_resolution_clock::now();
		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
		std::cout << std::endl << "Time: ~" << duration / 1000000 << " seconds." << std::endl << std::endl << std::endl;
	}

	return 0;
}