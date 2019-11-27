#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "sys/time.h"
#include "knnring.h"

#define MatrixOrder 1 // 1 = RowMajor, 2 = ColMajor, 3 = Check both

#define DefaultNumCorpusPoints 1423
#define DefaultDim 37
#define DefaultNumNeighbors 13


int n = DefaultNumCorpusPoints;
int d = DefaultDim;
int k = DefaultNumNeighbors;

static char *STR_CORRECT_WRONG[] = { (char*)"WRONG", (char*)"CORRECT" };
struct timeval startwtime, endwtime;

void help(int argc, char *argv[], int id);
double distColMajor(double *X, double *Y, int i, int j, int d, int n, int m);
double distRowMajor(double *X, double *Y, int i, int j, int d, int n, int m);
int validateResultColMajor(knnresult knnres, double * corpus, double * query, int n, int m, int d, int k);
int validateResultRowMajor(knnresult knnres, double * corpus, double * query, int n, int m, int d, int k);

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);       // initialize MPI
	int p, id;                    // # processess and PID
	MPI_Comm_rank(MPI_COMM_WORLD, &id); // Task ID
	MPI_Comm_size(MPI_COMM_WORLD, &p); // # tasks
	help(argc, argv, id);
	if (id == 0)               // ..... MASTER		
		printf("Running with values n=%i, d=%i, k=%i.\n", n, d, k);

	srand((unsigned int)time(NULL));
	double p_time;	// job time

	// Generate random point set
	double * corpusAll = NULL;
	double *corpus = NULL;
	if (id == 0) {                // ..... MASTER		
		printf("Generating random data set. ");
		gettimeofday(&startwtime, NULL);
		corpusAll = (double *)malloc(p*n*d * sizeof(double));
		for (int ip = 0; ip < p; ip++) {
			// "read" new chunk
			corpus = (double *)malloc(n*d * sizeof(double));
			for (int i = 0; i < n*d; i++) {
				corpusAll[i + ip * n*d] = ((double)(rand())) / (double)RAND_MAX;
				corpus[i] = corpusAll[i + ip * n*d];
			}
			if (ip == p - 1)            // last chunk is mine
				break;
			// which process to send? what tag?
			int dst = ip + 1;
			int tag = 1;
			// send to correct process
			MPI_Send(corpus, n*d, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
			free(corpus);
		}
		p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("DONE in %fsec!\n", p_time);
	}
	else {                      // ..... other processes
		// from which process to I receive (master)? what tag?
		int rcv = 0;
		int tag = 1;
		MPI_Status Stat;
		corpus = (double *)malloc(n*d * sizeof(double));
		MPI_Recv(corpus, n*d, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, &Stat);
	}

	// Find neighbors
	knnresult knnresall;
	knnresult knnres = distrAllkNN(corpus, n, d, k);
	// ~~~~~~~~~~~~~~~~~~~~ gather back kNN results
	if (id == 0) {                // ..... MASTER	
		printf("Finding all %i neighbors for all %i query points. ", k, n * p);
		gettimeofday(&startwtime, NULL);
		knnresall.nidx = (int *)malloc(n*p*k * sizeof(int));
		knnresall.ndist = (double *)malloc(n*p*k * sizeof(double));
		knnresall.m = n * p;
		knnresall.k = k;
		for (int ip = 0; ip < p - 1; ip++) {
			// from which process to I receive? what tag?
			int rcv = ip + 1;
			int tag = 1;
			MPI_Status Stat;
			MPI_Recv(&knnresall.nidx[ip*n*k], n*k, MPI_INT, rcv, tag, MPI_COMM_WORLD, &Stat);
			MPI_Recv(&knnresall.ndist[ip*n*k], n*k, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, &Stat);
		}
		// move my result to final struct
		for (int i = 0; i < n*k; i++) {
			knnresall.nidx[(p - 1)*n*k + i] = knnres.nidx[i];
			knnresall.ndist[(p - 1)*n*k + i] = knnres.ndist[i];
		}
		gettimeofday(&endwtime, NULL);
		p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("DONE in %fsec!\n", p_time);
	}
	else {                      // ..... other processes
		// which process to send? what tag?
		int dst = 0;
		int tag = 1;
		// send to correct process
		MPI_Send(knnres.nidx, n*k, MPI_INT, dst, tag, MPI_COMM_WORLD);
		MPI_Send(knnres.ndist, n*k, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
	}

	// Validate results
	if (id == 0) {                // ..... MASTER		
		printf("Validating results. ");
		gettimeofday(&startwtime, NULL);
		int isValid = (MatrixOrder == 1) ? validateResultRowMajor(knnresall, corpusAll, corpusAll, n*p, n*p, d, k)
			: (MatrixOrder == 2) ? validateResultColMajor(knnresall, corpusAll, corpusAll, n*p, n*p, d, k)
			: validateResultColMajor(knnresall, corpusAll, corpusAll, n*p, n*p, d, k) || validateResultRowMajor(knnresall, corpusAll, corpusAll, n*p, n*p, d, k);
		printf("Tester validation: %s NEIGHBORS. ", STR_CORRECT_WRONG[isValid]);
		gettimeofday(&endwtime, NULL);
		p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("DONE in %fsec!\n", p_time);
	}

	printf("Exiting from proccess id %i.\n", id);
	free(corpus);
	free(corpusAll);
	MPI_Finalize();               // clean-up
	return 0;
}

void help(int argc, char *argv[], int id)
{
	if (argc > 1) {
		for (int i = 1; i < argc; i += 2) {
			if (*argv[i] == '-') {
				if (*(argv[i] + 1) == 'n')
					n = atoi(argv[i + 1]);
				else if (*(argv[i] + 1) == 'd')
					d = atoi(argv[i + 1]);
				else if (*(argv[i] + 1) == 'k')
					k = atoi(argv[i + 1]);
				else {
					help(1, argv, id);
					return;
				}
			}
			else {
				help(1, argv, id);
				return;
			}
		}
		return;
	}
	if (id == 0) {
		printf("Flags to use:\n");
		printf("-n [Number] :Number of corpus points (default:%i)\n", DefaultNumCorpusPoints);
		printf("-d [Dimension] :Dimension of the space (default: %i)\n", DefaultDim);
		printf("-k [Number] :Number of neighbors (default:%i)\n", DefaultNumNeighbors);
	}
}

double distColMajor(double *X, double *Y, int i, int j, int d, int n, int m)
{
	/* compute distance */
	double dist = 0;
	for (int l = 0; l < d; l++) {
		dist += (X[l*n + i] - Y[l*m + j]) * (X[l*n + i] - Y[l*m + j]);
	}
	return sqrt(dist);
}

double distRowMajor(double *X, double *Y, int i, int j, int d, int n, int m)
{
	/* compute distance */
	double dist = 0;
	for (int l = 0; l < d; l++) {
		dist += (X[l + i * d] - Y[l + j * d]) * (X[l + i * d] - Y[l + j * d]);
	}
	return sqrt(dist);
}

int validateResultColMajor(knnresult knnres, double * corpus, double * query, int n, int m, int d, int k)
{
	/* loop through all query points */
	for (int j = 0; j < m; j++) {
		/* max distance so far (equal to kth neighbor after nested loop) */
		double maxDist = -1;
		/* mark all distances as not computed */
		int * visited = (int *)calloc(n, sizeof(int));
		/* go over reported k nearest neighbors */
		for (int i = 0; i < k; i++) {
			/* keep list of visited neighbors */
			visited[knnres.nidx[i*m + j]] = 1;
			/* get distance to stored index */
			double distxy = distColMajor(corpus, query, knnres.nidx[i*m + j], j, d, n, m);
			/* make sure reported distance is correct */
			if (fabs(knnres.ndist[i*m + j] - distxy) > 1e-8) return 0;
			/* distances should be non-decreasing */
			if (knnres.ndist[i*m + j] < maxDist) return 0;
			/* update max neighbor distance */
			maxDist = knnres.ndist[i*m + j];
		} /* for (k) -- reported nearest neighbors */
		/* now maxDist should have distance to kth neighbor */
		/* check all un-visited points */
		for (int i = 0; i < n; i++) {
			/* check only (n-k) non-visited nodes */
			if (!visited[i]) {
				/* get distance to unvisited vertex */
				double distxy = distColMajor(corpus, query, i, j, d, n, m);
				/* point cannot be closer than kth distance */
				if (distxy < maxDist) return 0;
			} /* if (!visited[i]) */
		} /* for (i) -- unvisited notes */
		/* deallocate memory */
		free(visited);
	} /* for (j) -- query points */
	return 1;
}

int validateResultRowMajor(knnresult knnres, double * corpus, double * query, int n, int m, int d, int k)
{
	/* loop through all query points */
	for (int j = 0; j < m; j++) {
		/* max distance so far (equal to kth neighbor after nested loop) */
		double maxDist = -1;
		/* mark all distances as not computed */
		int * visited = (int *)calloc(n, sizeof(int));
		/* go over reported k nearest neighbors */
		for (int i = 0; i < k; i++) {
			/* keep list of visited neighbors */
			visited[knnres.nidx[i + j * k]] = 1;
			/* get distance to stored index */
			double distxy = distRowMajor(corpus, query, knnres.nidx[i + j * k], j, d, n, m);
			/* make sure reported distance is correct */
			if (fabs(knnres.ndist[i + j * k] - distxy) > 1e-8) {
				printf("Error in dist: got %e instead of %e ", knnres.ndist[i + j * k], distxy);
				printf("at idy: %i with idx: %i. ", j, knnres.nidx[i + j * k]);
				return 0;
			}
			/* distances should be non-decreasing */
			if (knnres.ndist[i + j * k] < maxDist) {
				printf("Error in order: got %lf while max is %lf. ", knnres.ndist[i + j * k], maxDist);
				return 0;
			}
			/* update max neighbor distance */
			maxDist = knnres.ndist[i + j * k];
		} /* for (k) -- reported nearest neighbors */
		/* now maxDist should have distance to kth neighbor */
		/* check all un-visited points */
		for (int i = 0; i < n; i++) {
			/* check only (n-k) non-visited nodes */
			if (!visited[i]) {
				/* get distance to unvisited vertex */
				double distxy = distRowMajor(corpus, query, i, j, d, n, m);
				/* point cannot be closer than kth distance */
				if (distxy < maxDist) {
					printf("Error found other neighbor. ");
					return 0;
				}
			} /* if (!visited[i]) */
		} /* for (i) -- unvisited notes */
		/* deallocate memory */
		free(visited);
	} /* for (j) -- query points */
	return 1;
}
