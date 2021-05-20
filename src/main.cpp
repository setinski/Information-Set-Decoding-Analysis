#pragma omp parallel for
#pragma optimize( "", off )

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "space.h"
#include "isd.h"
#include "misc.h"

std::vector<unsigned int> ParseAlphabetSizes(std::string str) {
	std::stringstream ss = std::stringstream(str);
	std::vector<unsigned int> v;
	char ch;
	do {
		int i;
		ss >> i;
		v.push_back(i);
	} while (ss >> ch);
	return v;
}

int main()
{
	/* 1. Users' inputs */
	std::string metric, metricInput;
	do {
		std::cout << "METRIC CHOICE - Enter H for Hamming, ";
		std::cout << "L for Lee metric: ";
		std::cin >> metricInput;
		if (!metricInput.compare("H") || !metricInput.compare("h"))
			metric = "hamming";
		else if (!metricInput.compare("L") || !metricInput.compare("l"))
			metric = "lee";
		else
			metric = "";
	} while (!MetricCheck(metric));
	std::cout << std::endl;

	std::string algorithm, algorithmInput;
	do {
		std::cout << "ALGORITHM CHOICE - Enter P for Prange's algorithm, D ";
		std::cout << "for Dumer/Stern's algorithm, W for Wagner's algorithm: ";
		std::cin >> algorithmInput;
		if (!algorithmInput.compare("P") || !algorithmInput.compare("p"))
			algorithm = "prange";
		else if (!algorithmInput.compare("D") || !algorithmInput.compare("d"))
			algorithm = "dumer";
		else if (!algorithmInput.compare("W") || !algorithmInput.compare("w"))
			algorithm = "wagner";
		else
			algorithm = "";
	} while(!AlgCheck(algorithm));

	std::string quantumSwitch;
	do {
		std::cout << "Choose if the algorithms is quantum. Enter C ";
		std::cout << "if the algorithm is classical, Q if it is quantum: ";
		std::cin >> quantumSwitch;
	} while (quantumSwitch.compare("C") && quantumSwitch.compare("c")
			&& quantumSwitch.compare("Q") && quantumSwitch.compare("q"));
	if (!quantumSwitch.compare("C") || !quantumSwitch.compare("c"))
		quantum = 0;
	else
		quantum = 1;
	std::cout << std::endl;

	std::string sizes;
	std::cout << "ALPHABET SIZES - Choose the alphabet sizes.";
	std::cout << "Enter alphabet sizes separated by comas, ";
	std::cout << "for example: 3,5,7 (do not add space after comma): ";
	std::cin >> sizes;
	std::vector<unsigned int> alphabetSizes = ParseAlphabetSizes(sizes);
	std::cout << std::endl;

	std::cout << "SUMMARY: " << std::endl;
	std::cout << "metric: " << metric << std::endl;
	std::cout << "algorithm: " << algorithm << std::endl;
	std::cout << "quantum: " << quantum << std::endl;
	std::cout << std::endl;

	// 2. Calculating the hardest instances
	std::string outputFilename = "results";
	if (quantum)
	{
		outputFilename += "_quantum";
	}
	std::fstream outputFile(outputFilename + ".txt", std::ofstream::out |
													std::ofstream::trunc);
	outputFile << "alphabetSize codeRate weight optLevelNum paramL paramP";
	outputFile << "runtime(log 2) runtime(log alphabetSize)" << std::endl;
	outputFile.close();

	for (auto i = 0; i < alphabetSizes.size(); i++)
	{
		std::cout << "Prosessing alphabet size: " << alphabetSizes[i] << std::endl;
		auto t1 = std::chrono::high_resolution_clock::now();

		double paramL, paramP, weight, codeRate, runTime;
		unsigned int optLevelNum;
		std::function<double(double)> runTimeFun = [&](double cr) {return RunTime(metric, algorithm, alphabetSizes[i],
																			cr, weight, paramL, paramP, optLevelNum); };
		try
		{
			codeRate = GoldenSectionSearch(epsilon, 1.0 - epsilon, tol, runTimeFun);
		}
		catch (std::runtime_error& rte)
		{
			std::cout << rte.what() << std::endl;
			return -1;
		}
		catch (std::invalid_argument& ia)
		{
			std::cout << ia.what() << std::endl;
			return -1;
		}
		catch (...) { std::cout << "Default exception."; }
		
		try
		{
			runTime = RunTime(metric, algorithm, alphabetSizes[i], codeRate, weight, paramL, paramP, optLevelNum);
		}
		catch (std::runtime_error& rte)
		{
			std::cout << rte.what() << std::endl;
			return -1;
		}
		catch (std::invalid_argument& ia)
		{
			std::cout << ia.what() << std::endl;
			return -1;
		}
		catch (...) { std::cout << "Default exception."; }

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
		std::cout << std::endl << "Time elapsed: ~" << duration / 1000000 << " seconds." << std::endl << std::endl << std::endl;
	}

	std::cout << "FIINISHED - The parameters of the hardest instances can be found in ";
	std::cout << outputFilename << ".txt." << std::endl;

	return 0;
}