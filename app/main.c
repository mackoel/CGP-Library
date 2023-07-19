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

#define OPERATION_OPTB0 0
#define OPERATION_OPTB1 1
#define OPERATION_PRINT 2
#define OPERATION_TEST 3

char usage_string[256]= \
"Usage: cgpapp <mode> <parameters file>;\n" \ 
"cgpapp 1 parameters.txt train_filename_left train_filename_right\n" \
"cgpapp 3 parameters.txt test_filename_left test_filename_right 0 5\n";

int main(int argc, char** argv) {

	int operation_mode = 0; /* 0 - optimization
							   1 - print errors */

	if (argc > 1) {
		sscanf(argv[1], "%d", &operation_mode);
	} else {
		fprintf(stderr, "No parameters given, exiting\n");
		fprintf(stderr, usage_string);
		return 0;
	}

	struct parameters *paramsLeft = NULL;
	struct parameters *paramsRight = NULL;
	struct dataSet *trainingData = NULL;
	struct chromosome *chromo = NULL;

	int numInputsLeft = 4;
	int numInputsRight = 4;
	int numNodes = 200;
	int numOutputs = 1;
	int nodeArity = 2;
	int numGens = 80000;
	int numGensInt = 800;
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

	fscanf(pFile, "numInputsLeft = %i;\n", &numInputsLeft);
	fscanf(pFile, "numInputsRight = %i;\n", &numInputsRight);
	fscanf(pFile, "numNodes = %i;\n", &numNodes);
	fscanf(pFile, "numOutputs = %i;\n", &numOutputs);
	fscanf(pFile, "nodeArity = %i;\n", &nodeArity);
	fscanf(pFile, "numGens = %i;\n", &numGens);
	fscanf(pFile, "targetFitness = %lf;\n", &targetFitness);
	fscanf(pFile, "updateFrequency = %i;\n", &updateFrequency);
	fscanf(pFile, "maxMutationConst = %lf;\n", &maxMutationConst);
	fscanf(pFile, "numLevels = %i;\n", &numLevels);
	fscanf(pFile, "levelCoeff = %lf;\n", &levelCoeff);
	fscanf(pFile, "numGensInt = %i;\n", &numGensInt);

	double* defaultSimpleConstantsLeft = malloc(numInputsLeft * sizeof(double));

	fscanf(pFile, "\n", NULL);
	for (int i = 0; i < numInputsLeft; ++i) {
		fscanf(pFile, "%lf;\n", &defaultSimpleConstantsLeft[i]);
		// printf("%lf ", defaultSimpleConstants[i]);
	}

	fscanf(pFile, "\n", NULL);
	// printf("\n");
	double* shiftForSigmoidLeft = malloc(numInputsLeft * sizeof(double));
	for (int i = 0; i < numInputsLeft; ++i) {
		fscanf(pFile, "%lf;\n", &shiftForSigmoidLeft[i]);
		// printf("%lf ", shiftForSigmoid[i]);
	}

	fscanf(pFile, "\n", NULL);
	// printf("\n");
	double* scaleForSigmoidLeft = malloc(numInputsLeft * sizeof(double));
	for (int i = 0; i < numInputsLeft; ++i) {
		fscanf(pFile, "%lf;\n", &scaleForSigmoidLeft[i]);
		// printf("%lf ", scaleForSigmoid[i]);
	}

	double* defaultSimpleConstantsRight = malloc(numInputsRight * sizeof(double));

	for (int i = 0; i < numInputsRight; ++i) {
		fscanf(pFile, "%lf;\n", &defaultSimpleConstantsRight[i]);
		// printf("%lf ", defaultSimpleConstants[i]);
	}

	fscanf(pFile, "\n", NULL);
	// printf("\n");
	double* shiftForSigmoidRight = malloc(numInputsRight * sizeof(double));
	for (int i = 0; i < numInputsRight; ++i) {
		fscanf(pFile, "%lf;\n", &shiftForSigmoidRight[i]);
		// printf("%lf ", shiftForSigmoid[i]);
	}

	fscanf(pFile, "\n", NULL);
	// printf("\n");
	double* scaleForSigmoidRight = malloc(numInputsRight * sizeof(double));
	for (int i = 0; i < numInputsRight; ++i) {
		fscanf(pFile, "%lf;\n", &scaleForSigmoidRight[i]);
		// printf("%lf ", scaleForSigmoid[i]);
	}
	fclose(pFile);

	if (operation_mode == OPERATION_OPTB0) {

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
			return 0;
		}

		trainingData = initialiseDataSetFromFile(train_filename);

		paramsRight = initialiseParameters(numInputsRight, numNodes, numOutputs, nodeArity, maxMutationConst, defaultSimpleConstantsRight, shiftForSigmoidRight, scaleForSigmoidRight, -1);

		addNodeFunction(paramsRight, "add,sub,mul,div");
		addCustomNodeFunction(paramsRight, hyperbola, "hyperbola", 1);
		addCustomNodeFunction(paramsRight, linear, "linear", 1);

		setTargetFitness(paramsRight, targetFitness);

		setUpdateFrequency(paramsRight, updateFrequency);

		printParameters(paramsRight);

		for (int i = 0; i < numLevels; i++) {
			chromo = runCGP(paramsRight, trainingData, numGens);
		
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

			UpdateDataSet(paramsRight, chromo, trainingData, levelCoeff);

			freeChromosome(chromo);
		}

		freeParameters(paramsRight);
		freeDataSet(trainingData);
	}

	if (operation_mode == OPERATION_OPTB1) {

		char* train_filename_left = NULL;
		char* train_filename_right = NULL;
		char* test_filename = NULL;

		char* output_chromo = NULL;
		char* output_constants = NULL;
		char* output_tex = NULL;

		struct dataSet *trainingData_left = NULL;
		struct dataSet *trainingData_right = NULL;

		struct chromosome *chromo_left = NULL;
		struct chromosome *chromo_right = NULL;

		if (argc > 4) {
			train_filename_left = argv[3];
			train_filename_right = argv[4];
		} else {
			fprintf(stderr, "\n Incorrect input \n");
		}

		trainingData_left = initialiseDataSetFromFile(train_filename_left);
		trainingData_right = initialiseDataSetFromFile(train_filename_right);

		paramsLeft = initialiseParameters(numInputsLeft, numNodes, numOutputs, nodeArity, maxMutationConst, defaultSimpleConstantsLeft, shiftForSigmoidLeft, scaleForSigmoidLeft, -1);
		paramsRight = initialiseParameters(numInputsRight, numNodes, numOutputs, nodeArity, maxMutationConst, defaultSimpleConstantsRight, shiftForSigmoidRight, scaleForSigmoidRight, -1);
		
		addNodeFunction(paramsLeft, "add,sub,mul,div");
		addCustomNodeFunction(paramsLeft, hyperbola, "hyperbola", 1);
		addCustomNodeFunction(paramsLeft, linear, "linear", 1);

		addNodeFunction(paramsRight, "add,sub,mul,div");
		addCustomNodeFunction(paramsRight, hyperbola, "hyperbola", 1);
		addCustomNodeFunction(paramsRight, linear, "linear", 1);

		setTargetFitness(paramsLeft, targetFitness);
		setUpdateFrequency(paramsLeft, updateFrequency);
		printParameters(paramsLeft);

		setTargetFitness(paramsRight, targetFitness);
		setUpdateFrequency(paramsRight, updateFrequency);
		printParameters(paramsRight);

		chromo = runCGP(paramsRight, trainingData_right, numGens);
		
		// printChromosome(chromo, 0);

		char const_filename[100];
		char chromo_filename[100];
		char latex_filename[100];
		int i = 0;
		snprintf(const_filename, 100, "%s_const%02d.txt",train_filename_right, i);
		snprintf(chromo_filename, 100, "%s_chromo%02d.chromo", train_filename_right, i);
		snprintf(latex_filename, 100, "%s_latex%02d.tex", train_filename_right, i);
		saveConstants(chromo, const_filename);
		saveChromosome(chromo, chromo_filename);
		saveChromosomeLatex(chromo, 0, latex_filename);

		double* errors_chromo = (double*)calloc(getDataSetNumSamples(trainingData_right), sizeof(double));
		double* errors_chromo_left = (double*)calloc(getDataSetNumSamples(trainingData_left), sizeof(double));
		double* errors_chromo_right = (double*)calloc(getDataSetNumSamples(trainingData_right), sizeof(double));

		getResult(trainingData_right, errors_chromo, chromo, 1);

		updateDataSetOutput(trainingData_left, errors_chromo, levelCoeff);
		updateDataSetOutput(trainingData_right, errors_chromo, levelCoeff);

		// printChromosome(chromo, 0);

		freeChromosome(chromo);

		setUserData(paramsLeft, (void*)errors_chromo);
		setCustomFitnessFunction(paramsLeft, supervisedLearningUserData, "supervisedLearningUserData");

		setUserData(paramsRight, (void*)errors_chromo);
		setCustomFitnessFunction(paramsRight, supervisedLearningUserData, "supervisedLearningUserData");

		for (i = 1; i < numLevels; i++) {

			for (int j = 0; j < getDataSetNumSamples(trainingData_right); j++) {
				errors_chromo[j] = 1.0;
				errors_chromo_left[j] = 1.0;
				errors_chromo_right[j] = 1.0;
			}
			chromo_left = initialiseChromosome(paramsLeft);
			chromo_right = initialiseChromosome(paramsRight);

			int numIter = numGens / numGensInt;
			for (int k = 0; k < numIter; k++) {
				fprintf(stdout, "numLevel: %d of %d; ", i, numLevels);
				fprintf(stdout, "iter: %d of %d; left ", k, numIter);
				chromo_left = rerunCGP(paramsLeft, trainingData_left, numGensInt, chromo_left);
				getResult(trainingData_left, errors_chromo_left, chromo_left, 1);
				for (int j = 0; j < getDataSetNumSamples(trainingData_right); j++) {
					errors_chromo[j] = errors_chromo_left[j];
				}
				fprintf(stdout, "numLevel: %d of %d; ", i, numLevels);
				fprintf(stdout, "iter: %d of %d; right ", k, numIter);
				chromo_right = rerunCGP(paramsRight, trainingData_right, numGensInt, chromo_right);
				getResult(trainingData_right, errors_chromo_right, chromo_right, 1);
				for (int j = 0; j < getDataSetNumSamples(trainingData_right); j++) {
					errors_chromo[j] = errors_chromo_right[j];
				}
			}
			for (int j = 0; j < getDataSetNumSamples(trainingData_right); j++) {
				errors_chromo[j] = levelCoeff * errors_chromo_left[j] * errors_chromo_right[j];
			}
			updateDataSetOutput(trainingData_left, errors_chromo, 1);
			updateDataSetOutput(trainingData_right, errors_chromo, 1);
			
			// printChromosome(chromo_left, 0);

			snprintf(const_filename, 100, "%s_const_left%02d.txt",train_filename_left, i+1);
			snprintf(chromo_filename, 100, "%s_chromo_left%02d.chromo", train_filename_left, i+1);
			snprintf(latex_filename, 100, "%s_latex_left%02d.tex", train_filename_left, i+1);
			saveConstants(chromo_left, const_filename);
			saveChromosome(chromo_left, chromo_filename);
			saveChromosomeLatex(chromo_left, 0, latex_filename);

			// printChromosome(chromo_right, 0);

			snprintf(const_filename, 100, "%s_const_right%02d.txt",train_filename_right, i+1);
			snprintf(chromo_filename, 100, "%s_chromo_right%02d.chromo", train_filename_right, i+1);
			snprintf(latex_filename, 100, "%s_latex_right%02d.tex", train_filename_right, i+1);
			saveConstants(chromo_right, const_filename);
			saveChromosome(chromo_right, chromo_filename);
			saveChromosomeLatex(chromo_right, 0, latex_filename);

			freeChromosome(chromo_left);
			freeChromosome(chromo_right);
		}

		freeParameters(paramsLeft);
		freeParameters(paramsRight);
		freeDataSet(trainingData_left);
		freeDataSet(trainingData_right);
	}

	if (operation_mode == OPERATION_TEST) {

		char* train_filename_left = NULL;
		char* train_filename_right = NULL;
		char* test_filename = NULL;

		char* output_chromo = NULL;
		char* output_constants = NULL;
		char* output_tex = NULL;

		struct dataSet *trainingData_left = NULL;
		struct dataSet *trainingData_right = NULL;

		struct chromosome *chromo_left = NULL;
		struct chromosome *chromo_right = NULL;

		if (argc > 4) {
			train_filename_left = argv[3];
			train_filename_right = argv[4];
		} else {
			fprintf(stderr, "\n Incorrect input \n");
		}

		trainingData_left = initialiseDataSetFromFile(train_filename_left);
		trainingData_right = initialiseDataSetFromFile(train_filename_right);

		paramsLeft = initialiseParameters(numInputsLeft, numNodes, numOutputs, nodeArity, maxMutationConst, defaultSimpleConstantsLeft, shiftForSigmoidLeft, scaleForSigmoidLeft, -1);
		paramsRight = initialiseParameters(numInputsRight, numNodes, numOutputs, nodeArity, maxMutationConst, defaultSimpleConstantsRight, shiftForSigmoidRight, scaleForSigmoidRight, -1);

		addNodeFunction(paramsLeft, "add,sub,mul,div");
		addCustomNodeFunction(paramsLeft, hyperbola, "hyperbola", 1);
		addCustomNodeFunction(paramsLeft, linear, "linear", 1);

		addNodeFunction(paramsRight, "add,sub,mul,div");
		addCustomNodeFunction(paramsRight, hyperbola, "hyperbola", 1);
		addCustomNodeFunction(paramsRight, linear, "linear", 1);

		setTargetFitness(paramsLeft, targetFitness);
		setUpdateFrequency(paramsLeft, updateFrequency);
		// printParameters(paramsLeft);

		setTargetFitness(paramsRight, targetFitness);
		setUpdateFrequency(paramsRight, updateFrequency);
		// printParameters(paramsRight);

		char const_filename[100];
		char chromo_filename[100];
		char latex_filename[100];
		int i = 0;
		snprintf(const_filename, 100, "%s_const%02d.txt",train_filename_right, i);
		snprintf(chromo_filename, 100, "%s_chromo%02d.chromo", train_filename_right, i);
		snprintf(latex_filename, 100, "%s_latex%02d.tex", train_filename_right, i);

		chromo = initialiseChromosomeFromFileWithUserFunctions(paramsRight, chromo_filename);
		loadConstants(chromo, const_filename);
		
		// printChromosome(chromo, 0);

		double* errors_chromo = (double*)calloc(getDataSetNumSamples(trainingData_right), sizeof(double));
		double* errors_chromo_left = (double*)calloc(getDataSetNumSamples(trainingData_left), sizeof(double));
		double* errors_chromo_right = (double*)calloc(getDataSetNumSamples(trainingData_right), sizeof(double));

		getResult(trainingData_right, errors_chromo, chromo, levelCoeff);

		// updateDataSetOutput(trainingData_left, errors_chromo, 1);
		// updateDataSetOutput(trainingData_right, errors_chromo, 1);

		freeChromosome(chromo);

		for (i = 1; i < numLevels; i++) {

			snprintf(const_filename, 100, "%s_const_left%02d.txt",train_filename_left, i+1);
			snprintf(chromo_filename, 100, "%s_chromo_left%02d.chromo", train_filename_left, i+1);
			snprintf(latex_filename, 100, "%s_latex_left%02d.tex", train_filename_left, i+1);

			chromo_left = initialiseChromosomeFromFileWithUserFunctions(paramsLeft, chromo_filename);
			loadConstants(chromo_left, const_filename);

			// printChromosome(chromo_left, 0);

			getResult(trainingData_left, errors_chromo_left, chromo_left, 1);

			snprintf(const_filename, 100, "%s_const_right%02d.txt",train_filename_right, i+1);
			snprintf(chromo_filename, 100, "%s_chromo_right%02d.chromo", train_filename_right, i+1);
			snprintf(latex_filename, 100, "%s_latex_right%02d.tex", train_filename_right, i+1);

			chromo_right = initialiseChromosomeFromFileWithUserFunctions(paramsRight, chromo_filename);
			loadConstants(chromo_right, const_filename);

			// printChromosome(chromo_right, 0);

			getResult(trainingData_right, errors_chromo_right, chromo_right, 1);

			for (int j = 0; j < getDataSetNumSamples(trainingData_right); j++) {
				errors_chromo[j] += levelCoeff * errors_chromo_left[j] * errors_chromo_right[j];
			}

			freeChromosome(chromo_left);
			freeChromosome(chromo_right);
		}

		for (int j = 0; j < getDataSetNumSamples(trainingData_right); j++) {
			fprintf(stdout,"%4d, %13.6lf, %13.6lf\n", j, errors_chromo[j], getDataSetSampleOutput(trainingData_right, j, 0));
		}

		freeParameters(paramsLeft);
		freeParameters(paramsRight);
		freeDataSet(trainingData_left);
		freeDataSet(trainingData_right);
	}

	if (operation_mode == OPERATION_PRINT) {
		int begin = -1;
		int end = 18;

		char* test_filename = NULL;
		char* train_filename = NULL;
		if (argc == 7) {
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

			chromo = initialiseChromosomeFromFile(chromo_filename, maxMutationConst, defaultSimpleConstantsRight, shiftForSigmoidRight, scaleForSigmoidRight);
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