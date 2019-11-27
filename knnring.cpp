#include "pch.h"

#include "knnring.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cblas.h"
#include "mpi.h"


double *blas_distXY(double *X, double *Y, int n, int m, int d)
{
	// Calculate X.X**T
	double *dot_X = (double*)malloc(n * sizeof(double));
	for (int i = 0; i < n; i++)
		*(dot_X + i) = cblas_ddot(d, (X + i * d), 1, (X + i * d), 1);
	// Calculate Y.Y**T
	double *dot_Y = (double*)malloc(m * sizeof(double));
	for (int i = 0; i < m; i++)
		*(dot_Y + i) = cblas_ddot(d, (Y + i * d), 1, (Y + i * d), 1);
	// Calculate (X.X**T + Y.Y**T)
	double *D = (double*)calloc(m*n, sizeof(double));
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			*(D + i * n + j) = *(dot_Y + i) + *(dot_X + j);
	// Calculate 2.Y.X**T and all distances
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, d, -2.0, Y, d, X, d, 1.0, D, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++) {
			if (*(D + i * n + j) < 1e-12)	// threshold min value for openblas
				*(D + i * n + j) = 0.0;
			else
				*(D + i * n + j) = sqrt(*(D + i * n + j));
		}
	return D;
}

void SWAP(double *d, int *idx, int a, int b)
{
	double tmpd;
	tmpd = *(d + a);
	*(d + a) = *(d + b);
	*(d + b) = tmpd;
	int tmpi = *(idx + a);
	*(idx + a) = *(idx + b);
	*(idx + b) = tmpi;
}

double quick_select(double *d, int *idx, int len, int k)
{
	int i, st;
	for (st = i = 0; i < len - 1; i++) {
		if (d[i] > d[len - 1]) continue;
		SWAP(d, idx, i, st);
		st++;
	}
	SWAP(d, idx, len - 1, st);
	return k == st ? d[st]
		: st > k ? quick_select(d, idx, st, k)
		: quick_select(d + st, idx + st, len - st, k - st);
}

void my_qsort(double *d, int *idx, int len, int k)
{
	quick_select(d, idx, len, k);			// put first the k shortest
	for (int i = 0; i < k + 1; i++)			// and after short them
		quick_select(d, idx, k + 1, i);
}

void add_idx_offset(int offset, int *idx, int n, int k)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < k; j++)
			*(idx + i * k + j) += offset;
}

knnresult kNN(double *X, double  *Y, int n, int m, int d, int k)
{
	knnresult result;
	result.nidx = (int*)malloc(m * k * sizeof(int));
	result.ndist = (double*)malloc(m * k * sizeof(double));
	result.m = m;
	result.k = k;

	double *distxy = blas_distXY(X, Y, n, m, d);
	for (int q = 0; q < m; q++) {	// all query
		int *idx = (int*)malloc(n * sizeof(int));
		for (int c = 0; c < n; c++)	// create index array
			*(idx + c) = c;
		my_qsort((distxy + q * n), idx, n, k);
		for (int z = 0; z < k; z++) {	// add to struct the k neighbors
			*(result.nidx + q * k + z) = *(idx + z);
			*(result.ndist + q * k + z) = *(distxy + q * n + z);
		}
	}
	return result;
}

knnresult distrAllkNN(double *Y, int n, int d, int k)
{
	int p, id;							// # processess and PID
	MPI_Comm_rank(MPI_COMM_WORLD, &id); // Task ID
	MPI_Comm_size(MPI_COMM_WORLD, &p);	// # tasks
	// where receive? where send?
	int source = (id == 0) ? p - 1 : id - 1;	// if id == 0 take from last proccess
	int dest = (id == p - 1) ? 0 : id + 1;		// if last proccess send to first proccess
	int tag = (int)'X';					// mpi tag to send X data points set
	MPI_Status Stat;					// mpi receive status

	int *npidx = (int*)malloc(n * (2 * k) * sizeof(int));
	double *npdist = (double*)malloc(n * (2 * k) * sizeof(double));
	double *X = (double*)malloc(n * d * sizeof(double));
	double *X_recv = (double*)malloc(n * d * sizeof(double));
	for (int ip = 0; ip < p; ip++) {	// move points into circle
		if (ip == 0) {	// at first loop you already have data
			memcpy(X, Y, n * d * sizeof(double));
			knnresult new_knnresult;
			new_knnresult = kNN(X, Y, n, n, d, k);
			int offset = (id >= ip + 1) ? (id - ip - 1) * n : (p + id - ip - 1) * n;
			add_idx_offset(offset, new_knnresult.nidx, n, k);				// offset idx depending of the chunk
			for (int i = 0; i < n; i++) {
				memcpy((npidx + i * (2 * k)), (new_knnresult.nidx + i * k), k * sizeof(int));
				memcpy((npdist + i * (2 * k)), (new_knnresult.ndist + i * k), k * sizeof(double));
			}
		}
		else {	// work with received data
			knnresult new_knnresult;
			new_knnresult = kNN(X, Y, n, n, d, k);
			int offset = (id >= ip + 1) ? (id - ip - 1) * n : (p + id - ip - 1) * n;
			add_idx_offset(offset, new_knnresult.nidx, n, k);				// offset idx depending of the chunk
			for (int i = 0; i < n; i++) {
				memcpy((npidx + i * (2 * k) + k), (new_knnresult.nidx + i * k), k * sizeof(int));
				memcpy((npdist + i * (2 * k) + k), (new_knnresult.ndist + i * k), k * sizeof(double));
				my_qsort((npdist + i * (2 * k)), (npidx + i * (2 * k)), 2 * k, k);
			}
		}

		if (ip < p - 1) {	// don't send data at the last loop, the job was done at the first
			// send and receive data from other proccess in the ring
			if (id % 2) {	// even id number
				MPI_Send(X, n*d, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
				MPI_Recv(X_recv, n*d, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &Stat);
			}
			else {			// odd id number
				MPI_Recv(X_recv, n*d, MPI_DOUBLE, source, tag, MPI_COMM_WORLD, &Stat);
				MPI_Send(X, n*d, MPI_DOUBLE, dest, tag, MPI_COMM_WORLD);
			}
			memcpy(X, X_recv, n * d * sizeof(double));
		}
	}

	knnresult final_knnresult;
	final_knnresult.nidx = (int*)malloc(n * k * sizeof(int));
	final_knnresult.ndist = (double*)malloc(n * k * sizeof(double));
	for (int i = 0; i < n; i++) {
		memcpy((final_knnresult.nidx + i * k), (npidx + i * (2 * k)), k * sizeof(int));
		memcpy((final_knnresult.ndist + i * k), (npdist + i * (2 * k)), k * sizeof(double));
	}
	final_knnresult.m = n;
	final_knnresult.k = k;
	return final_knnresult;
}