#ifndef LIBCLUST_H_
#define LIBCLUST_H_

#include <stdlib.h>
#include <math.h>


double corr(double *obj1, double *obj2, int n);
double** corrDistMatrix(int nrows, int ncols, double** data);
double** compDistMatrix(int nrows, int ncols, double** data, int toPrint);

#endif

