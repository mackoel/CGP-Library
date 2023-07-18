#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgp.h"


double maximum(double num[], int size)
{
	int i;
	double max;
	max = num[0];

	for (i = 1; i < size; i++) {
		if (num[i] > max)
			max = num[i];
	}
	return max;
}


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


double hyperbola(const int numInputs, const double *inputs, const double *connectionWeights, double simpleConstant) {
	return 1 / (inputs[0] - simpleConstant);
}


double linear(const int numInputs, const double *inputs, const double *connectionWeights, double simpleConstant) {
	return (inputs[0] - simpleConstant);
}


double SS_tot(double a[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) sum += a[i];

	double mean = (double)sum /(double)n;

	double sqDiff = 0;
	for (int i = 0; i < n; i++) sqDiff += (a[i] - mean) * (a[i] - mean);
	return sqDiff;
}


double SS_res(double a[], int n)
{
	double sum = 0;
	for (int i = 0; i < n; i++) sum += a[i]*a[i];
	return sum;
}

int main(int argc, char** argv) {

	int operation_mode = 0; /* 0 - optimization
							   1 - print errors */

	if (argc > 1) {
		sscanf(argv[1], "%d", &operation_mode);
	} else {
		fprintf(stderr, "No parameters given, exiting\n");
		return 0;
	}

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
	int numLevels = 10;
	double levelCoeff = 0.9;

	FILE * pFile;
	pFile = fopen(argv[2], "r");
	if (pFile == NULL) {
		fprintf(stderr, "Can't open parameters file, exiting\n");
		return 0;
	}

	fscanf(pFile, "numInputs = %i;\n", &numInputs);
	fscanf(pFile, "numNodes = %i;\n", &numNodes);
	fscanf(pFile, "numOutputs = %i;\n", &numOutputs);
	fscanf(pFile, "nodeArity = %i;\n", &nodeArity);
	fscanf(pFile, "numGens = %i;\n", &numGens);
	fscanf(pFile, "targetFitness = %lf;\n", &targetFitness);
	fscanf(pFile, "updateFrequency = %i;\n", &updateFrequency);
	fscanf(pFile, "maxMutationConst = %lf;\n", &maxMutationConst);
	fscanf(pFile, "numLevels = %i;\n", &numLevels);
	fscanf(pFile, "levelCoeff = %lf;\n", &levelCoeff);

	double* defaultSimpleConstants = malloc(numInputs * sizeof(double));

	for (int i = 0; i < numInputs; ++i) {
		fscanf(pFile, "%lf;\n", &defaultSimpleConstants[i]);
		// printf("%lf ", defaultSimpleConstants[i]);
	}

	fscanf(pFile, "\n", NULL);
	printf("\n");
	double* shiftForSigmoid = malloc(numInputs * sizeof(double));
	for (int i = 0; i < numInputs; ++i) {
		fscanf(pFile, "%lf;\n", &shiftForSigmoid[i]);
		// printf("%lf ", shiftForSigmoid[i]);
	}

	fscanf(pFile, "\n", NULL);
	printf("\n");
	double* scaleForSigmoid = malloc(numInputs * sizeof(double));
	for (int i = 0; i < numInputs; ++i) {
		fscanf(pFile, "%lf;\n", &scaleForSigmoid[i]);
		// printf("%lf ", scaleForSigmoid[i]);
	}

	if (operation_mode == 0) {

		char* train_filename = NULL;
		char* test_filename = NULL;

		char* output_chromo = NULL;
		char* output_constants = NULL;
		char* output_tex = NULL;

		char* strNumOfSnp = NULL;
		if (argc > 3) {
			train_filename = argv[3];
		} else {
			fprintf(stderr, "\n Incorrect input \n");
		}

		trainingData = initialiseDataSetFromFile(train_filename);

		params = initialiseParameters(numInputs, numNodes, numOutputs, nodeArity, maxMutationConst, defaultSimpleConstants, shiftForSigmoid, scaleForSigmoid, -1);

		addNodeFunction(params, "add,sub,mul,div");
		addCustomNodeFunction(params, hyperbola, "hyperbola", 1);
		addCustomNodeFunction(params, linear, "linear", 1);

		setTargetFitness(params, targetFitness);

		setUpdateFrequency(params, updateFrequency);

		printParameters(params);

		for (int i = 0; i < numLevels; i++) {
			chromo = runCGP(params, trainingData, numGens);
		
			printChromosome(chromo, 0);

			char const_filename[100];
			char chromo_filename[100];
			char latex_filename[100];
			snprintf(const_filename, 100, "%s_const%02d.txt",train_filename, i+1); // i+1
			snprintf(chromo_filename, 100, "%s_chromo%02d.chromo", train_filename, i+1); //i+1
			snprintf(latex_filename, 100, "%s_latex%02d.tex", train_filename, i+1); //i+1
			saveConstants(chromo, const_filename);
			saveChromosome(chromo, chromo_filename);
			saveChromosomeLatex(chromo, 0, latex_filename);

			UpdateDataSet(params, chromo, trainingData, levelCoeff);

			freeChromosome(chromo);
		}

		freeParameters(params);
		freeDataSet(trainingData);
	}

	if (operation_mode == 1) {
		int begin = -1;
		int end = 18;

		char* test_filename = NULL;
		char* train_filename = NULL;
		if (argc == 5) {
			test_filename = argv[3];
			train_filename = argv[4];
			begin = -1;
			end = 18;
		} else if (argc == 7) {
			test_filename = argv[3];
			train_filename = argv[4];
			begin = strtol(argv[5], NULL, 10);
			end = strtol(argv[6], NULL, 10);
		} else {
			fprintf(stderr, "train.data test.data\n");
			return 0;
		}

		/*double error = 0;*/
		struct dataSet *data;
		data = initialiseDataSetFromFile(test_filename);
		double* errors = (double*)calloc(getDataSetNumSamples(data), sizeof(double));

		for (int i = begin; i < end; i++) {
			char const_filename[100];
			char chromo_filename[100];

			snprintf(const_filename, 100, "%s_const%02d.txt", train_filename, i+1); 
			snprintf(chromo_filename, 100, "%s_chromo%02d.chromo", train_filename, i+1); 

			chromo = initialiseChromosomeFromFile(chromo_filename, maxMutationConst, defaultSimpleConstants, shiftForSigmoid, scaleForSigmoid);
			loadConstants(chromo, const_filename);

			getResult(data, errors, chromo, levelCoeff);

			freeChromosome(chromo);
		}

		double* real = (double*)calloc(getDataSetNumSamples(data), sizeof(double));
		for (int i = 0; i < getDataSetNumSamples(data); i++)
			real[i] = getDataSetSampleOutput(data, i, 0);
		double SS_tot_ = SS_tot(real, getDataSetNumSamples(data));

		for (int i = 0; i < getDataSetNumSamples(data); i++)
			errors[i] = fabs(errors[i] - getDataSetSampleOutput(data, i, 0));

		double SS_res_ = SS_res(errors, getDataSetNumSamples(data));
		double sum_errors = 0;
		for (int i = 0; i < getDataSetNumSamples(data); i++)
			sum_errors += errors[i];

		printf("%lf\n", sum_errors);

		printf("%lf\n", 1 - SS_res_/SS_tot_);
		
		printf("%lf\n", maximum(errors, getDataSetNumSamples(data)));

		freeDataSet(data);
	}

	return 0;
}