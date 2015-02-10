#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "libClust.h"

double corr(double *obj1, double *obj2, int n)
/*
from: gsl_stats_correlation()
  Calculate Pearson correlation = cov(X, Y) / (sigma_X * sigma_Y)
This routine efficiently computes the correlation in one pass of the
data and makes use of the algorithm described in:

B. P. Welford, "Note on a Method for Calculating Corrected Sums of
Squares and Products", Technometrics, Vol 4, No 3, 1962.

This paper derives a numerically stable recurrence to compute a sum
of products

S = sum_{i=1..N} [ (x_i - mu_x) * (y_i - mu_y) ]

with the relation

S_n = S_{n-1} + ((n-1)/n) * (x_n - mu_x_{n-1}) * (y_n - mu_y_{n-1})
*/

{
    double sum_xsq;
    double sum_ysq;
    double sum_cross;
    double ratio;
    double delta_x, delta_y;
    double mean_x, mean_y;
    //double corrVal;

    int i;

    mean_x = obj1[0];
    mean_y = obj2[0];
    sum_xsq = 0.0;
    sum_ysq = 0.0;
    sum_cross = 0.0;

    /*
     * Compute:
     * sum_xsq = Sum [ (x_i - mu_x)^2 ],
     * sum_ysq = Sum [ (y_i - mu_y)^2 ] and
     * sum_cross = Sum [ (x_i - mu_x) * (y_i - mu_y) ]
     * using the above relation from Welford's paper
     */

    for(i = 1; i < n; i++)
    {
        ratio = i / (i + 1.0);
        delta_x = obj1[i] - mean_x;
        delta_y = obj2[i] - mean_y;
        sum_xsq += delta_x * delta_x * ratio;
        sum_ysq += delta_y * delta_y * ratio;
        sum_cross += delta_x * delta_y * ratio;
        mean_x += delta_x / (i + 1.0);
        mean_y += delta_y / (i + 1.0);
    }

    return sum_cross / (sqrt(sum_xsq) * sqrt(sum_ysq));
}

//==================================================================


double** corrDistMatrix(int nrows, int ncols, double** data)
/*
Purpose
=======

Calculate the distance matrix of a data matrix. It assumes that each row is an
object and the each column is a feature or variable. Furthermore it assumes that
there are no missing data.
The distance are based on the Pearson correlation formula.


Arguments
=========

nrows      (input) int
The number of rows of the data matrix

ncolumns   (input) int
The number of columns of the data matrix

data       (input) double[nrows][ncolumns]
The bidimensional array containing the data

Return
=========
matrix      (double**), distance matrix
*/
{
    int i,j;
    double** matrix;

    if (nrows < 2)
        return NULL;

    /* Set up the ragged array */
    matrix = malloc(nrows*sizeof(double*));
    if(matrix==NULL) return NULL; /* Not enough memory available */
    matrix[0] = NULL; /* The zeroth row has zero columns. We allocate it anyway for convenience.*/

    for (i = 1; i < nrows; i++)
    {
        matrix[i] = malloc(i*sizeof(double));
        if (matrix[i]==NULL)
            break; /* Not enough memory available */
    }
    if (i < nrows) /* break condition encountered */
    {
        j = i;
        for (i = 1; i < j; i++)
            free(matrix[i]);
        return NULL;
    }

    /* Calculate the distances and save them in the ragged array */
    for (i = 1; i < nrows; i++){
        j = 0;
        while (j < i){
            /* Compute distance */
            matrix[i][j] = corr(data[i], data[j], ncols);
            j++;
        }
     }
    return matrix;
}


double** compDistMatrix(int nrows, int ncols, double** data, int toPrint)
/*

Purpose
=======

Compute the distance matrix e print it if you want.
It's wrapper for corrDistMatrix.

Arguments
=========

nrows    int
number of rows of the dataset

ncols    int
number of cols of the dataset

data	double**
data matrix

toPrint int
enable (if it is 1) the check prints

Return
=========
distance matrix
*/
{
    int i, j;
    double** distMatrix;

    distMatrix = corrDistMatrix(nrows, ncols, data);
    if (!distMatrix)
    {
        printf ("Insufficient memory to store the distance matrix\n");
        return NULL;
    }

    if (toPrint == 1){
        printf("\n============ Distance matrix ============\n");
        printf("   Obj:");
        for(i=0; i<nrows-1; i++)
            printf("%6d", i);
        printf("\n");
        for(i=0; i<nrows; i++)
        {
            printf("Obj %2d:",i);
            for(j=0; j<i; j++) printf(" %5.2f",distMatrix[i][j]);
                printf("\n");
        }
        printf("\n");
    }
    return distMatrix;
}
