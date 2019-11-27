#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "sys/time.h"
#include "knnring.h"

#define MatrixOrder 1 // 1 = RowMajor, 2 = ColMajor, 3 = Check both

#define DefaultNumCorpusPoints 897
#define DefaultNumQueryPoints 762
#define DefaultDim 7
#define DefaultNumNeighbors 13


int n = DefaultNumCorpusPoints;
int m = DefaultNumQueryPoints;
int d = DefaultDim;
int k = DefaultNumNeighbors;

static char *STR_CORRECT_WRONG[] = { (char*)"WRONG", (char*)"CORRECT" };
struct timeval startwtime, endwtime;

void help(int argc, char *argv[]);
double distColMajor(double *X, double *Y, int i, int j, int d, int n, int m);
double distRowMajor(double *X, double *Y, int i, int j, int d, int n, int m);
int validateResultColMajor(knnresult knnres, double * corpus, double * query, int n, int m, int d, int k);
int validateResultRowMajor(knnresult knnres, double * corpus, double * query, int n, int m, int d, int k);

int main(int argc, char *argv[])
{
	help(argc, argv);
	printf("Running with values n=%i, m=%i, d=%i, k=%i.\n", n, m, d, k);

	srand((unsigned int)time(NULL));

	// Generate random point set
	printf("Generating random data set. ");
	gettimeofday(&startwtime, NULL);
	double *corpus = (double *)malloc(n * d * sizeof(double));
	double *query = (double *)malloc(m * d * sizeof(double));
	for (int i = 0; i < n*d; i++)
		*(corpus + i) = ((double)rand() / (RAND_MAX));
	for (int i = 0; i < m*d; i++)
		*(query + i) = ((double)rand() / (RAND_MAX));
	gettimeofday(&endwtime, NULL);
	double p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("DONE in %fsec!\n", p_time);

	// Find neighbors
	printf("Finding all %i neighbors for all %i query points. ", k, m);
	gettimeofday(&startwtime, NULL);
	knnresult knnres = kNN(corpus, query, n, m, d, k);
	gettimeofday(&endwtime, NULL);
	p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("DONE in %fsec!\n", p_time);

	// Validate results
	printf("Validating results. ");
	gettimeofday(&startwtime, NULL);
	int isValid = (MatrixOrder == 1) ? validateResultRowMajor(knnres, corpus, query, n, m, d, k)
		: (MatrixOrder == 2) ? validateResultColMajor(knnres, corpus, query, n, m, d, k)
		: validateResultColMajor(knnres, corpus, query, n, m, d, k) || validateResultRowMajor(knnres, corpus, query, n, m, d, k);
	printf("Tester validation: %s NEIGHBORS. ", STR_CORRECT_WRONG[isValid]);
	gettimeofday(&endwtime, NULL);
	p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("DONE in %fsec!\n", p_time);

	printf("Exiting\n");
	free(query);
	free(corpus);
	return 0;
}


void help(int argc, char *argv[])
{
	if (argc > 1) {
		for (int i = 1; i < argc; i += 2) {
			if (*argv[i] == '-') {
				if (*(argv[i] + 1) == 'n')
					n = atoi(argv[i + 1]);
				else if (*(argv[i] + 1) == 'm')
					m = atoi(argv[i + 1]);
				else if (*(argv[i] + 1) == 'd')
					d = atoi(argv[i + 1]);
				else if (*(argv[i] + 1) == 'k')
					k = atoi(argv[i + 1]);
				else {
					help(1, argv);
					return;
				}
			}
			else {
				help(1, argv);
				return;
			}
		}
		return;
	}
	printf("Flags to use:\n");
	printf("-n [Number] :Number of corpus points (default:%i)\n", DefaultNumCorpusPoints);
	printf("-m [Number] :Number of query points (default:%i)\n", DefaultNumQueryPoints);
	printf("-k [Number] :Number of neighbors (default:%i)\n", DefaultNumNeighbors);
	printf("-d [Dimension] :Dimension of the space (default: %i)\n", DefaultDim);
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
