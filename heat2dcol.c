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
	double *x;
	double **w;
	double **full;
	double *ghost1;
	double *ghost2;
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
		printf("HEAT2D Col Decomposition\n");
		printf("  C version\n");
		printf("  A program to solve for the steady state temperature distribution\n");
		printf("  over a rectangular plate using col decomposition.\n");
		printf("  Spatial grid of %d by %d points.\n", M, N);
		printf("\n");
	}


	// Calculate array size for each process
	if (my_id == 0) {
		u = calcArraySize(M, (N / nproc) + (N % nproc) + 1);
		w = calcArraySize(M, (N / nproc) + (N % nproc) + 1);
		full = calcArraySize(M, N);
	} else if (my_id == (nproc - 1)) {
		u = calcArraySize(M, (N / nproc) + 1);
		w = calcArraySize(M, (N / nproc) + 1);
	}  else {
		u = calcArraySize(M, (N / nproc) + 2);
		w = calcArraySize(M, (N / nproc) + 2);
	}


	// Used to convert column vector to row vector
	ghost1 = malloc(M * sizeof(double *));
	ghost2 = malloc(M * sizeof(double *));	
	x = malloc(M * (N / nproc) * sizeof(double *));


	// Print if process 0
	if (my_id == 0) {
		printf("\n");
		printf("  The iteration will be repeated until the change is <= %G\n", eps);
		printf("  Boundary Temperatures  left: %G  right: %G  top: %G  bottom: %G\n", Tl, Tr, Tt, Tb);
		printf("  The steady state solution will be written to '%s'.\n", output_file);
	}


	diff = eps;

	
 	// Set the left/right boundary values, which don't change.
	if (my_id == 0 || my_id == (nproc - 1)) {
		for (i = 0; i < M; i++) {
			if (my_id == 0) {
				w[i][0] = Tl;
			} else {
				w[i][(N / nproc)] = Tr;
			}
		}
	}


	// Set the top/bottom boundary values, which don't change.
	if (my_id == 0) {
		for (j = 0; j < ((N / nproc) + (N % nproc) + 1); j++) {
			w[0][j] = Tt;
			w[M - 1][j] = Tb;
		}
	} else if (my_id == (nproc - 1)) {
		for (j = 0; j < ((N / nproc) + 1); j++) {
			w[0][j] = Tt;
			w[M - 1][j] = Tb;
		}
	} else {
		for (j = 0; j < ((N / nproc) + 2); j++) {
			w[0][j] = Tt;
			w[M - 1][j] = Tb;
		}
	}
	
	
 	// Average the boundary values, to come up with a reasonable initial value for the interior
	mean = ((N * Tt) + (N * Tb) + ((M - 2) * Tl) + ((M - 2) * Tr)) / (2 * M + 2 * N - 4);


 	// Initialize the interior solution to the mean value.
 	if (my_id == 0) {
		w = fillInteriorCells(1, 1, M - 1, (N / nproc) + (N % nproc) + 1, mean, w);
	} else if (my_id == (nproc - 1)){
		w = fillInteriorCells(1, 0, M - 1, (N / nproc), mean, w);
	} else {
		w = fillInteriorCells(1, 0, M - 1, (N / nproc), mean, w);
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
			u = copyArray(M, (N / nproc) + (N % nproc) + 1, w, u);
		} else if (my_id == (nproc - 1)) {
			u = copyArray(M, (N / nproc) + 1, w, u);
		} else {
			u = copyArray(M, (N / nproc) + 2, w, u);
		}
	
	
		// Determine the new estimate of the solution at the interior points.
		// The new solution W is the average of north, south, east and west neighobrs.
		diff = 0.0;
		if (my_id == 0) {
			diff = calcDiff(M - 1, (N / nproc) + (N % nproc), w, u, diff);
		} else if (my_id == (nproc - 1)) {
			diff = calcDiff(M - 1, (N / nproc), w, u, diff);
		} else {
			diff = calcDiff(M - 1, (N / nproc) + 1, w, u, diff);
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


		// Pass common columns to neighbor processes
		if (nproc > 1) {
			if (eps <= maxDiff) {
				if (my_id == 0) {
					for (i = 0; i < M; i++) {
						ghost1[i] = w[i][(N / nproc) + (N % nproc) - 1];
					}
					MPI_Isend(ghost1, M, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
				} else if (my_id == (nproc - 1)) {
					for (i = 0; i < M; i++) {
						ghost1[i] = w[i][1];
					}
					MPI_Isend(ghost1, M, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
				} else {
					for (i = 0; i < M; i++) {
						ghost1[i] = w[i][1];
						ghost2[i] = w[i][(N / nproc)];
					}
					MPI_Isend(ghost1, M, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
					MPI_Isend(ghost2, M, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
				}

				MPI_Barrier(MPI_COMM_WORLD);
		
				if (my_id == 0) {
                			MPI_Irecv(ghost1, M, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
					for (i = 0; i < M; i++) {
						w[i][(N / nproc) + (N % nproc)] = ghost1[i];
					}
                		} else if (my_id == (nproc - 1)) {
                        		MPI_Irecv(ghost1, M, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
                			for (i = 0; i < M; i++) {
						w[i][0] = ghost1[i];
					}
				} else {
                        		MPI_Irecv(ghost1, M, MPI_DOUBLE, my_id - 1, 0, MPI_COMM_WORLD, &request);
                        		MPI_Irecv(ghost2, M, MPI_DOUBLE, my_id + 1, 0, MPI_COMM_WORLD, &request);
                			for (i = 0; i < M; i++) {
						w[i][0] = ghost1[i];
						w[i][(N / nproc) + 1] = ghost2[i];
					}
				}
				
				// Block all processes until the receives have complete
				MPI_Barrier(MPI_COMM_WORLD);	
			}
		}
	} while (eps <= maxDiff);


	// Finish clock
	ctime2 = cpu_time();
	ctime = ctime2 - ctime1;
	

	// Send process data to process 0
	if (my_id != 0) {

		// Transfer 2D array to a 1D array
		int a = 0;
		for (i = 0; i < M; i++) {
			for (j = 1; j < (N / nproc) + 1; j++) {
				x[a] = w[i][j];
				a++;
			}
		}

		MPI_Send(x, (N / nproc) * M, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}


	// Gather process data and store to file
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
		
		// Transfer process 0 to M x N array
		for (i = 0; i < M; i++) {
			for (j = 0; j < ((N / nproc) + (N % nproc)); j++) {
				full[i][j] = w[i][j];
			}
		}


		// Transfer process 1 to (nproc - 1) to M x N array
		int p;
		int k = (N / nproc) + (N % nproc);
		for (p = 1; p < nproc; p++) {
			MPI_Recv(x, (N / nproc) * M, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &status);
			int l = 0;
			for (i = 0; i < M; i++) {
				for (j = k; j < k + (N / nproc); j++) {
					full[i][j] = x[l];
					l++;
				}
			}
			k += (N / nproc);
		}

		
		// Write M x N array to file
		for (i = 0; i < M; i++) {
			for (j = 0; j < N; j++) {
				fprintf(fp, "%6.2f ", full[i][j]);
			}
			fputc('\n', fp);
		}

		fclose(fp);

		printf("\n");
		printf("  Solution written to the output file '%s'\n", output_file);

		// All done
		printf("\n");
		printf("HEAT2DCOL:\n");
		printf("  Normal end of execution.\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
