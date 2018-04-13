#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "svalidate.h"
#include "heat2dhelper.h"

int main(int argc, char *argv[]) {

	int i;
	int j;
	int M;
	int N;
	int iterations;
	int iterations_print;
	int success;
	int validData;
	int nproc;
	int my_id;
	float eps;
	double Tl, Tr, Tt, Tb;
	double diff;
	double maxDiff;
	double mean;
	double **u;
	double **x;
	double **w;
	double ctime;
	double ctime1;
	double ctime2;
	char *output_file;
	FILE *fp;
	MPI_Status status;
	MPI_Request request;


	// Initialize MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	
	// Verify command line input is correct format
	validData = verifyInput(argc, argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
	if (validData == 0) {
		if (my_id == 0) {
			usageRow();
		}
		MPI_Finalize();
		exit(-1);
	}


	// Initialize command line input
	M = atoi(argv[1]);
	N = atoi(argv[2]);
	Tl = atof(argv[3]);
	Tr = atof(argv[4]);
	Tt = atof(argv[5]);
	Tb = atof(argv[6]);
	eps = atof(argv[7]);
        output_file = argv[8];
	

	// Print if process 0
	if (my_id == 0) {
		printf("HEAT2D Row Decomposition\n");
		printf("  C version\n");
		printf("  A program to solve for the steady state temperature distribution\n");
		printf("  over a rectangular plate using row decomposition.\n");
		printf("  Spatial grid of %d by %d points.\n", M, N);
		printf("\n");
	}


	// Calculate array size for each process
	if (my_id == 0) {
		u = calcArraySize((M / nproc) + (M % nproc) + 1, N);
		w = calcArraySize((M / nproc) + (M % nproc) + 1, N);
	} else if (my_id == (nproc - 1)) {
		u = calcArraySize((M / nproc) + 1, N);
		w = calcArraySize((M / nproc) + 1, N);
	}  else {
		u = calcArraySize((M / nproc) + 2, N);
		w = calcArraySize((M / nproc) + 2, N);
	}


	// Print if process 0
	if (my_id == 0) {
		printf("\n");
		printf("  The iteration will be repeated until the change is <= %G\n", eps);
		printf("  Boundary Temperatures  left: %G  right: %G  top: %G  bottom: %G\n", Tl, Tr, Tt, Tb);
		printf("  The steady state solution will be written to '%s'.\n", output_file);
	}


	diff = eps;

	
 	// Set the left/right boundary values, which don't change.
	if (my_id == 0) {
		for (i = 0; i < ((M / nproc) + (M % nproc) + 1); i++) {
			w[i][0] = Tl;
			w[i][N-1] = Tr;
		}
	} else if (my_id == (nproc - 1)) {
		for (i = 0; i < ((M / nproc) + 1); i++) {
			w[i][0] = Tl;
			w[i][N - 1] = Tr;
		}
	} else {
		for (i = 0; i < ((M / nproc) + 2); i++) {
			w[i][0] = Tl;
			w[i][N - 1] = Tr;
		}
	}

	// Set the top/bottom boundary values, which don't change.
	if (my_id == 0 || my_id == (nproc - 1)) {
		for (j = 0; j < N; j++) {
			if (my_id == 0) {
				w[0][j] = Tt;
			} else {
				w[(M / nproc)][j] = Tb;
			}
		}
	}

	
 	// Average the boundary values, to come up with a reasonable initial value for the interior
	mean = ((N * Tt) + (N * Tb) + ((M - 2) * Tl) + ((M - 2) * Tr)) / (2 * M + 2 * N - 4);


 	// Initialize the interior solution to the mean value.
 	if (my_id == 0) {
		w = fillInteriorCells(1, 1, (M / nproc) + (M % nproc) + 1, N - 1, mean, w);
	} else if (my_id == (nproc - 1)){
		w = fillInteriorCells(0, 1, (M / nproc), N - 1, mean, w);
	} else {
		w = fillInteriorCells(0, 1, (M / nproc) + 2, N - 1, mean, w);
	}

	
 	// Iterate until the new solution W differs form the old solution U by no more than EPSILON.
 	iterations = 0;
	iterations_print = 1;

	if (my_id == 0) {
		printf("\n");
		printf("  Iteration  Change\n");
		printf("\n");
	}

	
	// Start clock
	ctime1 = cpu_time();


	// Main calculation
	do {
		// Save the old solution in u
		if (my_id == 0) {
			u = copyArray((M / nproc) + (M % nproc) + 1, N, w, u);
		} else if (my_id == (nproc - 1)) {
			u = copyArray((M / nproc) + 1, N, w, u);
		} else {
			u = copyArray((M / nproc) + 2, N, w, u);
		}
		
		
		// Determine the new estimate of the solution at the interior points.
		// The new solution W is the average of north, south, east and west neighobrs.
		diff = 0.0;
		if (my_id == 0) {
			diff = calcDiff((M / nproc) + (M % nproc), N - 1, w, u, diff);
		} else if (my_id == (nproc - 1)) {
			diff = calcDiff((M / nproc), N - 1, w, u, diff);
		} else {
			diff = calcDiff((M / nproc) + 1, N - 1, w, u, diff);
		}

		
		// Block all processes until averaging is done
		MPI_Barrier(MPI_COMM_WORLD);

		// Get max diff from all processes
		MPI_Allreduce(&diff, &maxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		
		// Block all processes until maxDiff is calculated (probably delete)
		MPI_Barrier(MPI_COMM_WORLD);

		
		iterations++;
		if (iterations == iterations_print) {
			if (my_id == 0) {
				printf("  %8d  %f\n", iterations, maxDiff);
			}

			iterations_print = 2 * iterations_print;
		}

		
		// Pass common rows to neighbor processes
		if (nproc > 1) {
			if (eps <= maxDiff) {
				if (my_id == 0) {
					MPI_Isend(w[(M / nproc) + (M % nproc) - 1], N, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
				} else if (my_id == (nproc - 1)) {
					MPI_Isend(w[1], N, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
				} else {
					MPI_Isend(w[1], N, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
					MPI_Isend(w[(M / nproc)], N, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
				}

				MPI_Barrier(MPI_COMM_WORLD);
	
				if (my_id == 0) {
                			MPI_Irecv(w[(M / nproc) + (M % nproc)], N, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
                		} else if (my_id == (nproc - 1)) {
                        		MPI_Irecv(w[0], N, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
                		} else {
                        		MPI_Irecv(w[0], N, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
                        		MPI_Irecv(w[(M / nproc) + 1], N, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
                		}	
			}
		}

		// Block all processes until passing common rows completes
		MPI_Barrier(MPI_COMM_WORLD);

	} while (eps <= maxDiff);


	// Finish clock
	ctime2 = cpu_time();
	ctime = ctime2 - ctime1;
	

	// Send process data to process 0
	if (my_id != 0) {
		MPI_Send(&w[1][0], 2 * (M / nproc) * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}


	// Gather data in process 0 and save to file
	if (my_id == 0) {
		printf("\n");
		printf("  %8d  %f\n", iterations, maxDiff);
		printf("\n");
		printf("  Error tolerance achieved.\n");
		printf("  CPU time = %f\n", ctime);


		// Write the solution to the output file
		fp = fopen(output_file, "w");

		fprintf(fp, "%d\n", M);
		fprintf(fp, "%d\n", N);
	
		// Write process 0 to file
		for (i = 0; i < (M / nproc) + (M % nproc); i++) {
			for (j = 0; j < N; j++) {
				fprintf(fp, "%6.2f ", w[i][j]);
			}
			fputc('\n', fp);
		}
	
		// Write process 1 to (nproc - 1) to file
		int p;
		for (p = 1; p < nproc; p++) {
			MPI_Recv(&u[0][0], 2 * (M / nproc) * N, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);
			for (i = 0; i < (M / nproc); i++) {
				for (j = 0; j < N; j++) {
					fprintf(fp, "%6.2f ", u[i][j]);
				}
				fputc('\n', fp);
			}
		}

		fclose(fp);

		printf("\n");
		printf("  Solution written to the output file '%s'\n", output_file);

		// All done
		printf("\n");
		printf("HEAT2DROW:\n");
		printf("  Normal end of execution.\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
