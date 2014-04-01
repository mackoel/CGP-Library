/*
	This file is part of CGP-Library
	Copyright (c) Andrew James Turner 2014 (andrew.turner@york.ac.uk)

    CGP-Library is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published 
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CGP-Library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with CGP-Library.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/cgp.h"  


float meanSquareError(struct parameters *params, struct chromosome *chromo, struct dataSet *data){
	
	int i,j;
	float squareError = 0;
			
	/* error checking */
	if(getNumChromosomeInputs(chromo) !=getNumDataSetInputs(data)){
		printf("Error: the number of chromosome inputs must match the number of inputs specified in the dataSet.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	if(getNumChromosomeOutputs(chromo) != getNumDataSetOutputs(data)){
		printf("Error: the number of chromosome outputs must match the number of outputs specified in the dataSet.\n");
		printf("Terminating CGP-Library.\n");
		exit(0);
	}

	/* for each sample in data */
	for(i=0; i<getNumDataSetSamples(data); i++){
	
		/* calculate the chromosome outputs for the set of inputs  */
		executeChromosome(chromo, getDataSetSampleInputs(data, i));
	
		/* for each chromosome output */
		for(j=0; j<getNumChromosomeOutputs(chromo); j++){
				
			squareError += powf(getDataSetSampleOutput(data,i,j) - getChromosomeOutput(chromo,j), 2);
		}
	}
	
	return squareError / (getNumDataSetSamples(data) * getNumDataSetOutputs(data));
}


int main(void){
	
	struct parameters *params;
	
	int numInputs = 2;
	int numNodes = 10;
	int numOutputs = 1;
	int arity = 3;
		
	params = initialiseParameters(numInputs, numNodes, numOutputs, arity);
	
	setRandomNumberSeed(123456789);	
	
	addNodeFunction(params, "add,sub,mul,sq,sin");

	setFitnessFunction(params, meanSquareError, "MSE");
	
	printParameters(params);
	
	struct dataSet *data = initialiseDataSetFromFile("./examples/symbolic.data");
	
	struct chromosome *chromo = runCGP(params, data, 10000);
	
	
	freeChromosome(chromo);
	freeDataSet(data);
	freeParameters(params);
	
	return 1;
	
	return 1;
}









