#include <stdio.h>
#include <stdlib.h>  
#include <time.h>
#include "cluster.h" /* The C Clustering Library */
#include "libwwd.h"


/* ========================================================================= */

double** distanceMatrix(int nrows, int ncols, double** data, int** mask, char metric, int toPrint)
/*

Purpose
=======

Compute the distance matrix

Arguments
=========

nrows    int
number of rows of the dataset

ncols    int
number of cols of the dataset

data	double**
data matrix

mask	double**
matrix that cointains information about missing data

metric	char
metric or distance to use for computing the distance matrix

toPrint int
enable (if it is 1) the check prints

========================================================================
*/
{ 
	int i, j;
	double** distMatrix;
	double* weight = malloc(ncols*sizeof(double));

    for (i = 0; i < ncols; i++)
		weight[i] = 1.0;
	distMatrix = distancematrix(nrows, ncols, data, mask, weight, metric, 0);

	if (!distMatrix)
	{ 
		printf ("Insufficient memory to store the distance matrix\n");
		free(weight);
		return NULL;
	}
	if (toPrint == 1){
        printf("\n============ Distance Matrix ============\n");
        printf("   Obj:");
        for(i=0; i<nrows-1; i++)
            printf("%6d", i);
        printf("\n");
        for(i=0; i<nrows; i++)
        {
            printf("Obj %2d:",i);
            for(j=0; j<i; j++)
                printf(" %5.2f",distMatrix[i][j]);
            printf("\n");
        }
     }
	
	free(weight);
	return distMatrix;
}


/* ========================================================================= */

void hclust(int nrows, int ncols, double** distmatrix, int toPrint, int numClusters)
/*

Purpose
=======

Compute the hierarchical clustering

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

numClusters int
number of clusters to cut
========================================================================
*/
{
    int i;
    const int nnodes = nrows-1;
    double* weight = malloc(ncols*sizeof(double));
    int* clusterid;
    Node* tree;
	
    if (numClusters > nrows){
        printf("The number of clusters is greater than the number of rows\n");
        return;
    }
	
    for (i = 0; i < ncols; i++) weight[i] = 1.0;

    /* Hierachical clustering, default distance is euclidean */
    tree = palcluster(nrows, distmatrix);

    //tree = treecluster(nrows, ncols, data, mask, weight, 0, metric, method, distmatrix);
    if (!tree)
    {
        /* Indication that the treecluster routine failed */
        printf ("treecluster routine failed due to insufficient memory\n");
        free(weight);
        return;
    }
/*	printf("Node     Item 1   Item 2    Distance\n");*/
/*	for(i=0; i<nnodes; i++)*/
/*	printf("%3d:%9d%9d      %g\n",*/
/*		   -i-1, tree[i].left, tree[i].right, tree[i].distance);*/
/*	printf("\n");*/

    /* Cutting the tree */
    if (numClusters != -1){
        clusterid = malloc(nrows*sizeof(int));
        cuttree (nrows, tree, numClusters, clusterid);


    if (toPrint == 1){
        for(i=0; i<nrows; i++)
            printf("Gene %2d: cluster %2d\n", i, clusterid[i]);
        printf("\n");
    }
        free(clusterid);
    }
    free(tree);
    free(weight);
    return;
}


Node *hAvgClust(int nrows, double** distmatrix, int toPrint, int numClusters)
/*

Purpose
=======

Compute the hierarchical clustering with average method

Arguments
=========

nrows    int
number of rows of the dataset

data	distmatrix**
distance matrix

toPrint int
enable (if it is 1) the check prints

numClusters int
number of clusters for cutting (currently it don't do cutting, add cuttree function if you want it)


Return
======
tree    Node*
Hierachical tree
========================================================================
*/
{
    Node* tree;
    int k;

    if (numClusters > nrows){
        printf("The number of clusters is greater than the number of rows\n");
        return;
    }


    /* Hierachical clustering, default distance is euclidean */
    tree = palcluster(nrows, distmatrix);

    if (!tree)
    {
        /* Indication that the treecluster routine failed */
        printf ("treecluster routine failed due to insufficient memory\n");
        return;
    }

    if (toPrint == 1){
        printf("\n==== hierarchical clustering output ==== \n");
        printf("Node     Item 1   Item 2    Distance\n");
        for(k=0; k<nrows-1; k++)
            printf("%3d:%9d%9d      %g\n", -k-1, tree[k].left, tree[k].right, tree[k].distance);
        printf("\n");
    }

    //cuttree here if you want
    return tree;
}



Node* hclustering(int nrows, int ncols, double** data, double** distmatrix, char distance, char method, int toPrint)
/*
 * Purpose
 * =======
 * Compute hierachical clustering. It's a wrapper for treecluster();
 *
 * Arguments
 * =========
 * nrows: (int), number of rows of data. It' also the number of objects;
 * ncols: (int), number of columns of data.
 * data: (double**)[nrows, ncols], matrix of data.
 * distmatrix: (double**)[nrows, nrows], matrix of distances of data.
 * distance: (char) type of distance to use for calculating distmatrix if it is not passed.
 * method: (char), method of linkage to use in the clustering
 * toPrint: (int),if 1 print clustering output
 *
 * return
 * ======
 * tree: (Node*), output of hierarchical clustering.
 */
{
    int i,j;
    double* weight = malloc(ncols*sizeof(double));
    int** mask = malloc(nrows*sizeof(double*));
    Node* tree;

    for (i = 0; i < ncols; i++){
        weight[i] = 1.0;
    }

    for (i=0;i<nrows;i++)
        mask[i] = malloc(ncols*sizeof(int));

    for (i=0;i<nrows;i++)
        for (j=0;j<ncols; j++)
            mask[i][j]=1;

    tree = treecluster(nrows, ncols, data, mask, weight, 0, distance, method, distmatrix);
    if (!tree)
    { /* Indication that the treecluster routine failed */
        printf ("treecluster routine failed due to insufficient memory\n");
        free(weight);
        return;
    }

    if (toPrint == 1){
        printf("\n==== hierarchical clustering output ==== \n");
        printf("Node     Item 1   Item 2    Distance\n");
        for(i=0; i<nrows-1; i++)
            printf("%3d:%9d%9d      %g\n", -i-1, tree[i].left, tree[i].right, tree[i].distance);
        printf("\n");
    }

    for (i = 0; i < nrows; i++){
        free(mask[i]);
    }
    free(mask);
    free(weight);

    return tree;
}




Node* hclusteringM(int nrows, int ncols, double** data, int** mask, double** distmatrix, char distance, char method, int toPrint)
/*
 * Purpose
 * =======
 * Compute hierachical clustering. It's a wrapper for treecluster();
 *
 * Arguments
 * =========
 * nrows: (int), number of rows of data. It' also the number of objects;
 * ncols: (int), number of columns of data.
 * data: (double**)[nrows, ncols], matrix of data.
 * mask: (int**), missing data mask of data
 * distmatrix: (double**)[nrows, nrows], matrix of distances of data.
 * distance: (char) type of distance to use for calculating distmatrix if it is not passed.
 * method: (char), method of linkage to use in the clustering
 * toPrint: (int),if 1 print clustering output
 *
 * return
 * ======
 * tree: (Node*), output of hierarchical clustering.
 */
{
    int i,j;
    double* weight = malloc(ncols*sizeof(double));
    Node* tree;

    for (i = 0; i < ncols; i++){
        weight[i] = 1.0;
    }

    tree = treecluster(nrows, ncols, data, mask, weight, 0, distance, method, distmatrix);
    if (!tree)
    { /* Indication that the treecluster routine failed */
        printf ("treecluster routine failed due to insufficient memory\n");
        free(weight);
        return;
    }

    if (toPrint == 1){
        printf("\n==== hierarchical clustering output ==== \n");
        printf("Node     Item 1   Item 2    Distance\n");
        for(i=0; i<nrows-1; i++)
            printf("%3d:%9d%9d      %g\n", -i-1, tree[i].left, tree[i].right, tree[i].distance);
        printf("\n");
    }

    free(weight);

    return tree;
}
