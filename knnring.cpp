#include "pch.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "knnring.h"


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
	quick_select(d, idx, len, k);
	for (int i = 0; i < k + 1; i++)
		quick_select(d, idx, k + 1, i);
}

knnresult kNN(double *X, double  *Y, int n, int m, int d, int k) 
{
	knnresult result;
	result.nidx = (int*)malloc(m * k * sizeof(int));
	result.ndist = (double*)malloc(m * k * sizeof(double));
	result.m = m;
	result.k = k;

	for (int q = 0; q < m; q++) {	// all query
		double *distxy = (double*)malloc(n * sizeof(double));
		int *idx = (int*)malloc(n * sizeof(int));
		for (int c = 0; c < n; c++) {	// all corpus
			*(idx + c) = c;
			*(distxy + c) = 0.0;
			for (int dim = 0; dim < d; dim++)
				*(distxy + c) += pow(*(Y + q * d + dim) - *(X + c * d + dim), 2);
			*(distxy + c) = sqrt(*(distxy + c));			
		}
		my_qsort(distxy, idx, n, k);
		for (int n = 0; n < k; n++) {	// add to struct the k neighbors
			*(result.nidx + n) = *(idx + n);
			*(result.ndist + n) = *(distxy + n);
		}
	}
	return result;
}

knnresult distrAllkNN(double *X, int n, int d, int k)
{
	knnresult result;
	result.nidx = (int*)malloc(n*d * sizeof(int));
	result.ndist = (double*)malloc(n*d * sizeof(double));
	result.m = n;
	result.k = k;

	return result;
}