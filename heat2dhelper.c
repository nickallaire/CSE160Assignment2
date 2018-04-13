#include "heat2dhelper.h"
#include "svalidate.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>


/*
 * \fn cpu_time
 * \param none
 * \return calculates the current reading of the CPU clock and returns in double format
 */
double cpu_time(void) {
	double value;
	value = (double) clock() / (double) CLOCKS_PER_SEC;
	return value;
}


/*
 * \fn usageRow()
 * \param none
 * \prints usage error if incorrect input 
 */
int usageRow() {
	fprintf(stderr, "usage: heat2drow M N Tl Tr Tt Tb eps file\n");
}


/*
 * \fn usageCol()
 * \param none
 * \prints usage error is incorrect input
 */
int usageCol() {
	fprintf(stderr, "usage: heat2dcol M N Tl Tr Tt Tb eps file\n");
}


/*
 * \fn verifyInput()
 * \param argc - number of command line arguments, c1 - M input, c2 - N input, c3 - Tl input, 
 *        c4 - Tr input, c5 - Tt input, c6 - Tb input, c7 - eps input
 * \brief checks user input from command line and verifies everything is in 
 *        the proper format
 * \return 1 if all data is valid, 0 otherwise
 */
int verifyInput(int argc, char c1[], char c2[], char c3[], char c4[], char c5[], char c6[], char c7[]) {
	if (argc < 9) {
		//printf("Not enough arguments\n");
                return 0;
        }
	
        if (!isInteger(trim(c1)) || atoi(trim(c1)) < 10) {
		//printf("Can't convert M to int or incorrect value\n");
                return 0;
        }

        if (!isInteger(trim(c2)) || atoi(trim(c2)) < 10) {
                //printf("Can't convert N to int or incorrect value\n");
		return 0;
        }

        if (!isDouble(trim(c3))) {
		//printf("Can't convert Tl to double\n");
                return 0;
        }

        if (!isDouble(trim(c4))) {
		//printf("Can't convert Tr to double\n");
                return 0;
        }

        if (!isDouble(trim(c5))) {
		//printf("Can't convert Tt to double\n");
                return 0;
        }

        if (!isDouble(trim(c6))) {
		//printf("Can't convert Tb to double\n");
                return 0;
        }

	if (!isFloat(trim(c7)) || atof(trim(c7)) >= 1.0 || atof(trim(c7)) < 0.0) {
		//printf("Can't convert eps to float or incorrect value\n");
                return 0;
        }

	return 1;
}


/*
 * \fn calcArraySize
 * \param row - number of rows to allocate, col - number of columns to allocate
 * \brief mallocs proper rows and columns to create a 2D array of pointers
 * \return array of size row x col
 */
double **calcArraySize(int row, int col) {
	double **array;
	int i;
	
	// Allocate rows
	array = (double **) malloc(row * sizeof(double *));
	
	// Allocate columns
	for (i = 0; i < row; i++) {
		array[i] = (double *) malloc(col * sizeof(double));
	}

	return array;
}


/*
 * \fn fillInteriorCells
 * \input startI - starting index for outter for loop, startJ - starting index for inner for loop, 
 *        row - number of rows to loop through, col - number of columns to loop through, 
 *        mean - mean of border values, array - array to fill in with mean value
 * \return an array with the interior cells filled with the mean
 */
double **fillInteriorCells(int startI, int startJ, int row, int col, double mean, double **array) {
	int i, j;

	for (i = startI; i < row; i++) {
		for (j = startJ; j < col; j++) {
			array[i][j] = mean;
		}
	}

	return array;
}


/*
 * \fn copyArray
 * \param row - number of rows in array, col - number of columns in array, 
 *        arrayFrom - array to copy from, arrayTo - array to copy to
 * \return array with copied values
 */
double **copyArray(int row, int col, double **arrayFrom, double **arrayTo) {
	int i, j;

	for (i = 0; i < row; i++) {
		for (j = 0; j < col; j++) {
			arrayTo[i][j] = arrayFrom[i][j];
		}
	}

	return arrayTo;
}


/*
 * \fn calcDiff
 * \param row - number of rows in array, col - number of columns in array,
 *        array - array to change values in, arrayOld - array with old values, 
 *        diff - absolute value difference between cell in array and arrayOld
 * \brief calculate the average value of bordering cells and store in array. 
 *        calculate the difference between new average and old average and 
 *        check if the absolute average of the difference if greater than current
 *        difference.
 * \return the maximum difference of a cell change in all cells in array
 */
double calcDiff(int row, int col, double **array, double **arrayOld, double diff) {
	int i, j;
	double oldVal;

	for (i = 1; i < row; i++) {
		for (j = 1; j < col; j++) {
			oldVal = array[i][j];
			array[i][j] = (arrayOld[i-1][j] + arrayOld[i+1][j] + arrayOld[i][j-1] + arrayOld[i][j+1]) / 4.0;
			if (diff < fabs(array[i][j] - arrayOld[i][j])) {
				diff = fabs(array[i][j] - arrayOld[i][j]);
			} 
		}
	}

	return diff;
}
