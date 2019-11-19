#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sys/time.h"
#include "knnring.h"

#define DefaultNumCorpusPoints 897
#define DefaultNumQueryPoints 762
#define DefaultDim 7
#define DefaultNumNeighbors 13


int n = DefaultNumCorpusPoints;
int m = DefaultNumQueryPoints;
int d = DefaultDim;
int k = DefaultNumNeighbors;
_Bool matlab = 0;

struct timeval startwtime, endwtime;

void help(int argc, char *argv[]);
void export_data(FILE *file, double *X, int x, int y, char *name);
void export_struct(FILE *file, knnresult knnres);

int main(int argc, char *argv[])
{
	help(argc, argv);
	printf("Running with values n=%i, m=%i, d=%i, k=%i.\n", n, m, d, k);

	srand((unsigned int)time(NULL));
	FILE *data = NULL;

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

	if (matlab) {
		printf("Writting data set to data.m. ");
		gettimeofday(&startwtime, NULL);
		data = fopen("./data.m", "w");
		fprintf(data, "n = %i;\n", n);
		fprintf(data, "m = %i;\n", m);
		fprintf(data, "d = %i;\n", d);
		fprintf(data, "k = %i;\n", k);
		export_data(data, corpus, n, d, (char*)"X");
		export_data(data, query, m, d, (char*)"Y");
		fclose(data);
		gettimeofday(&endwtime, NULL);
		p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("DONE in %fsec!\n", p_time);
	}

	// Find neighbors
	printf("Finding all %i neighbors for all %i query points. ", k, m);
	gettimeofday(&startwtime, NULL);
	knnresult knnres = kNN(corpus, query, n, m, d, k);
	gettimeofday(&endwtime, NULL);
	p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
	printf("DONE in %fsec!\n", p_time);

	if (matlab) {
		printf("Appending tree to data.m ");
		gettimeofday(&startwtime, NULL);
		data = fopen("./data.m", "a");
		export_struct(data, knnres);
		fclose(data);
		gettimeofday(&endwtime, NULL);
		p_time = (double)((endwtime.tv_usec - startwtime.tv_usec) / 1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
		printf("DONE in %fsec!\n", p_time);
	}

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
				else if (*(argv[i] + 1) == 'o') {
					matlab = 1;
					i--;
				}
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
	printf("-o :Print results into data.m file to evaluate in MATLAB\n");
}

void export_data(FILE *file, double *X, int x, int y, char *name)
{
	fprintf(file, "%s=[", name);
	for (int i = 0; i < x; i++) {
		fprintf(file, "[");
		for (int j = 0; j < y; j++)
			fprintf(file, "%lf ", *(X + i * y + j));
		fprintf(file, "]; ");
	}
	fprintf(file, "];\n");
}

void export_struct(FILE *file, knnresult knnres)
{
	for (int i = 0; i < knnres.m; i++) {
		fprintf(file, "result.list%i=[", i + 1);
		for (int j = 0; j < knnres.k; j++)
		{
			fprintf(file, "[%i %lf]; ", *(knnres.nidx + i * knnres.k + j) + 1, *(knnres.ndist + i * knnres.k + j));
		}
		fprintf(file, "];\n");
	}
}
