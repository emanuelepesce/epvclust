#ifndef DRIVERCLUSTER_H_
#define DRIVERCLUSTER_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cluster.h" /* The C Clustering Library */
#include "libwwd.h"

double** distanceMatrix(int nrows, int ncols, double** data, int** mask, char metric, int toPrint);
void hclust(int nrows, int ncols, double** distmatrix, int toPrint, int numClusters);
Node *hAvgClust(int nrows, double** distmatrix, int toPrint, int numClusters);
Node* hclustering(int nrows, int ncols, double** data, double** distmatrix, char method, char distance, int toPrint);
Node* hclusteringM(int nrows, int ncols, double** data, int** mask, double** distmatrix, char distance, char method, int toPrint);

#endif


