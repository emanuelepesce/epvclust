#ifndef PVCLUST_H_
#define PVCLUST_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include "pvclust.h"
#include "libwwd.h"
#include "libsort.h"
#include "driverCluster.h"
#include "libClust.h"
#include <limits.h>

void pvclustInit(char* filename, int nboot, int *r, int lenR);
int genRandInt(int limSup);
double** buildSample(char* filename, int dim);
double** getMinor(double **distMat, int dim, int* indices);
void freeMat(double **mat, int nRows);
void freeMatInt(int **mat, int nRows);
int** eboot(int toPrint);
void printVariables();
void epvclust(char* filename, int nboot);
int** buildPattern();

#endif

