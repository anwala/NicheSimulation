/*
	Niche Simulation:
	Simulation of mutation accumulation in the spermatogonial stem cell niche. 
	Code simulates the occurrence of mutations and the positive selection of 
	mutant cells within the niche over non-mutant cells.
	Stochastic progression of cell division is averaged over a number of 
	niches and outputs to file.

	Eoin C Whelan*, Alexander C Nwala, Christopher Osgood, Stephan Olariu
	Old Dominion University Norfolk VA USA 
*/
#include <iostream>
#include <time.h>
#include <random>

#include <sstream>
#include <fstream>
using namespace std;

std::default_random_engine generator(time(NULL));

struct Parameters
{
	int n;
	int divisions;
	int mutants;
	int runs;

	long double p;
	long double r;

	int *runsArray;
	float *averageArray;
	
};

/*
	Responsible for writing an experiment report to file

	@param arrayOfItems: list of things to be written to file
	@param lengthOfArray: size of arrayOfItems
	@param schema: descriptive information about what is being written
	@param filename: name of output file
*/
void writingToFileFunction(float* arrayOfItems, int lengthOfArray, string schema, string filename = "filename")
{
	ofstream myfile;
	filename = filename + ".csv";
	myfile.open(filename);

	myfile << "\"" << schema.c_str() << "\"" << endl;
	for (int i = 0; i < lengthOfArray; i++)
	{
		myfile << arrayOfItems[i] << endl;
	}

	myfile.close();

}

/*
	Generates a real random number between 0 and 1

	@return: a random number between 0 and 1
*/
long double rand01()
{
	std::uniform_real_distribution<long double> distribution(0.0, 1.0);
	long double randomBetween0and1 = distribution(generator);

	return randomBetween0and1;
}

/*
	Responsible running a simulation based on parameters

	@param experimentParameters: experiment parameters (e.g p, r)
	@return: the number of mutants
*/
int nicheSimulation(Parameters experimentParameters)
{
	// cout << "START" << endl;
	if (experimentParameters.n > 0 && experimentParameters.divisions > 0 && experimentParameters.p >= 0 && experimentParameters.r >= 0)
	{

		long double mutantProportion = experimentParameters.mutants / long double(experimentParameters.n);
		//cout << "mutants:" << mutants << endl;

		long double randomNumberBetween0and1 = rand01();
		//1 cell division
		if (randomNumberBetween0and1 < mutantProportion)
		{
			//cout << "...mutant divide" << endl;
			//this means we've picked a single mutant
			randomNumberBetween0and1 = rand01();

			if (randomNumberBetween0and1 < experimentParameters.r)
			{
				//cout << "......positive selection" << endl;
				//positive selection 
				//extract a random cell
				randomNumberBetween0and1 = rand01();
				mutantProportion = (experimentParameters.mutants + 1) / long double(experimentParameters.n + 1);
				if (randomNumberBetween0and1 < mutantProportion)
				{
					//cout << ".........mutant cell ejected" << endl;
				}
				else
				{
					//cout << ".........wild type cell ejected" << endl;
					//a wild type cell has been ejected
					experimentParameters.mutants++;
				}
			}
		}
		else
		{
			//cout << "...wild type cell divides" << endl;
			//we've picked a wild type cell
			randomNumberBetween0and1 = rand01();
			if (randomNumberBetween0and1 < experimentParameters.p)
			{
				//cout << "......mutant cell number incremented" << endl;

				experimentParameters.mutants++;
			}
		}
		
	}

	return experimentParameters.mutants;
}

/*
	Responsible for running multiple experiments

	@param: experiment parameters (e.g p, r)
	@return: an array in which every entry consists of the average number of mutant for a single cell division over all the runs
*/
float* runNExperiments(Parameters experimentParameters)
{
	if (experimentParameters.runs > 0 && experimentParameters.divisions > 0)
	{
		int mutantCount = 0;
		int i = 0;
		int j = 0;

		for (i = 0; i < experimentParameters.divisions; i++)
		{
			cout << "Division " << i << " of " << experimentParameters.divisions << endl;

			int totalMutantCellCount = 0;
			//this is for a single cell division instance run multiple independent times
			for (j = 0; j < experimentParameters.runs; j++)
			{
				experimentParameters.mutants = experimentParameters.runsArray[j];
				mutantCount = nicheSimulation(experimentParameters);
				experimentParameters.runsArray[j] = mutantCount;

				totalMutantCellCount += mutantCount;
			}

			//average mutant count for a single cell division over all the runs
			//in the event that experimentParameters.averageArray is too large to reside in memory,
			//it can be written directly to file at this point
			experimentParameters.averageArray[i] = totalMutantCellCount / float(experimentParameters.runs);
			//cout << "cell division:" << i + 1 << endl;
		}	
	}

	return experimentParameters.averageArray;
}

int main()
{
	Parameters experimentParameters;

	experimentParameters.n = 50;
	experimentParameters.divisions = 100000;
	experimentParameters.runs = 100;
	experimentParameters.p = 0.00005;
	experimentParameters.r = 0.01;
	experimentParameters.mutants = 0;

	experimentParameters.runsArray = new int[experimentParameters.runs];
	for (int i = 0; i < experimentParameters.runs; i++)
	{
		experimentParameters.runsArray[i] = 0;
	}

	experimentParameters.averageArray = new float[experimentParameters.divisions];
	for (int i = 0; i < experimentParameters.divisions; i++)
	{
		experimentParameters.averageArray[i] = 0;
	}


	//initialize
	clock_t tStart = clock();

	string schema =
		"n: " + to_string(experimentParameters.n) + "\n" +
		"divisions: " + to_string(experimentParameters.divisions) + "\n" +
		"runs: " + to_string(experimentParameters.runs) + "\n" +
		"p: " + to_string(experimentParameters.p) + "\n" +
		"r: " + to_string(experimentParameters.r) + "\n" +
		": " + to_string(experimentParameters.mutants) + "\n\n" +
		"AverageMutantCount";

	float *averageArray = runNExperiments(experimentParameters);

	writingToFileFunction(averageArray, experimentParameters.divisions, schema, "bigUrn_test");

	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	cout << "DONE" << endl;








	return 0;
}		