#include <stdio.h>
#include <stdlib.h>
#include "cgp.h"



int cmp_ptr(const void *a, const void *b)
{
	const int **left = (const int **)a;
	const int **right = (const int **)b;

	return (**left < **right) - (**right < **left);
}

size_t * order_int(const int *a, size_t n)
{
	const int **pointers = malloc(n * sizeof(const int *));
	for (size_t i = 0; i < n; i++) pointers[i] = a + i;
	qsort(pointers, n, sizeof(const int *), cmp_ptr);
	size_t *indices = malloc(n * sizeof(size_t));
	for (size_t i = 0; i < n; i++) indices[i] = pointers[i] - a;
	free(pointers);
	return indices;
}


double my_func(const int numInputs, const double *inputs, const double *connectionWeights, double simpleConstant) {
	return 1 / (inputs[0] - simpleConstant);
}

double my_func1(const int numInputs, const double *inputs, const double *connectionWeights, double simpleConstant) {
	return (inputs[0] - simpleConstant);
}






int main(int argc, char** argv) {


	struct parameters *params = NULL;
	struct dataSet *trainingData = NULL;
	struct chromosome *chromo = NULL;


	int numInputs = 4;
	int numNodes = 200;
	int numOutputs = 1;
	int nodeArity = 2;
	int numGens = 80000;
	double targetFitness = 0.01;
	int updateFrequency = 500;
	double maxMutationConst = 0.1;

	FILE * pFile;
	pFile = fopen("parameters.txt", "r");

	fscanf(pFile, "numInputs = %i;\n", &numInputs);
	fscanf(pFile, "numNodes = %i;\n", &numNodes);
	fscanf(pFile, "numOutputs = %i;\n", &numOutputs);
	fscanf(pFile, "nodeArity = %i;\n", &nodeArity);
	fscanf(pFile, "numGens = %i;\n", &numGens);
	fscanf(pFile, "targetFitness = %lf;\n", &targetFitness);
	fscanf(pFile, "updateFrequency = %i;\n", &updateFrequency);
	fscanf(pFile, "maxMutationConst = %lf;\n", &maxMutationConst);


	double* defaultSimpleConstants = malloc(numInputs * sizeof(double));

	for (int i = 0; i < numInputs; ++i) {
		fscanf(pFile, "%lf;\n", &defaultSimpleConstants[i]);
		printf("%lf ", defaultSimpleConstants[i]);
	}
	//defaultSimpleConstants[0] = 50000;
	//defaultSimpleConstants[1] = 0;
	//defaultSimpleConstants[2] = 0;
	//defaultSimpleConstants[3] = 0;
	fscanf(pFile, "\n", NULL);
	printf("\n");
	double* shiftForSigmoid = malloc(numInputs * sizeof(double));
	for (int i = 0; i < numInputs; ++i) {
		fscanf(pFile, "%lf;\n", &shiftForSigmoid[i]);
		printf("%lf ", shiftForSigmoid[i]);
	}

	fscanf(pFile, "\n", NULL);
	printf("\n");
	double* scaleForSigmoid = malloc(numInputs * sizeof(double));
	for (int i = 0; i < numInputs; ++i) {
		fscanf(pFile, "%lf;\n", &scaleForSigmoid[i]);
		printf("%lf ", scaleForSigmoid[i]);
	}


	// Note: you may need to check this path such that it is relative to your executable 
	//trainingData = initialiseDataSetFromFile("hardExample.data");

	char* train_filename = NULL;
	char* test_filename = NULL;

	char* output_chromo = NULL;
	char* output_constants = NULL;
	char* output_tex = NULL;

	char* strNumOfSnp = NULL;
	if (argc == 2) {
		train_filename = argv[1];
	}
	else {
		printf("\n Incorrect input \n");
	}

	trainingData = initialiseDataSetFromFile(train_filename);


	//int numOfSNP = atoi(strNumOfSnp);



		params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity, maxMutationConst, defaultSimpleConstants, shiftForSigmoid, scaleForSigmoid, -1);

		addNodeFunction(params, "add,sub,mul,div");
		addCustomNodeFunction(params, my_func, "hyperbola", 1);
		addCustomNodeFunction(params, my_func1, "linear", 1);

		setTargetFitness(params, targetFitness);

		setUpdateFrequency(params, updateFrequency);

		printParameters(params);

		chromo = runCGP(params, trainingData, numGens);
	
		printChromosome(chromo, 0);

		char const_filename[100];
		char chromo_filename[100];
		char latex_filename[100];
		snprintf(const_filename, 100, "%s_const%02d.txt",train_filename, 0); // i+1
		snprintf(chromo_filename, 100, "%s_chromo%02d.chromo", train_filename, 0); //i+1
		snprintf(latex_filename, 100, "%s_latex%02d.tex", train_filename, 0); //i+1
		saveConstants(chromo, const_filename);
		saveChromosome(chromo, chromo_filename);
		saveChromosomeLatex(chromo, 0, latex_filename);

		UpdateDataSet(params, chromo, trainingData);

		freeChromosome(chromo);
		freeParameters(params);




	freeDataSet(trainingData);


	return 0;
}