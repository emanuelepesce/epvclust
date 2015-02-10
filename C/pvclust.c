/*
* epvclust main
*
* Optimization in C of pvlcust algorithm
*
* Author: Emanuele Pesce
* Tutors: prof. Roberto Tagliaferri, PhD Angela Serra
*/
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

/*
 * Global variables
 */
static int nBoot;           // number of bootstrap replications
static int *dimsBoot;       // specifies the relative sample sizes of bootstrap replications
static int lenDimsBoot;     // lenght of dimsBoot
static char *filen;         // name of the dataset file
static int rowData;         // number of rows of the dataset
static int colData;         // number of columns of the dataset



void pvclustInit(char* filename, int nboot, int *r, int lenR)
/*
 * Purpose
 * =======
 * Inizializes variables;
 *
 * Arguments
 * =========
 * filename: (char*), name of dataset file;
 * nboot: (int), number of bootstrap replications;
 * r: (int*)[1,lenR], numeric vector that specifies the relative sample sizes of bootstrap replications;
 * lenR: (int), length of r;
 */
{
    int *dim,i;
    FILE *filep;

    filep = fopen(filename, "r");
    if (filep != NULL){
        /*	Get the number of columns and rows of the dataset	*/
        dim = getDatasetSize(filename);
        colData = dim[0];   //It assumes that objects are on the columns
        rowData = dim[1];
        free(dim);
        fclose(filep);
    }
    else{
        perror(filename);
        return;
    }

    filen = malloc(strlen(filename) + 1);
    strcpy(filen, filename);

    if (nboot > 0)
        nBoot = nboot;
    else
        nBoot = 1000;

    if ((r != 0) && (lenR>0)){
        dimsBoot = r;
        lenDimsBoot = lenR;
    }
    else{
        lenDimsBoot = 5;
        dimsBoot = malloc(sizeof(int)*lenDimsBoot);
        dimsBoot[0] = 50;
        for (i = 1; i < lenDimsBoot; i++){
            dimsBoot[i] = dimsBoot[i-1]+25;
        }
    }
}



void printVariables()
/*
 * Purpose
 * =======
 * Print variables. Useful for testing;
 *
 */
{
    int i;

    printf("\n===== Variables =====");
    printf("\nFilename: %s", filen);
    printf("\nnboot: %d", nBoot);
    printf("\nrowData: %d", rowData);
    printf("\ncolData: %d", colData);
    printf("\nlenDimsBoot: %d", lenDimsBoot);
    printf("\ndimsBoot :");
    for (i = 0; i < lenDimsBoot; i++){
        printf(" %d ", dimsBoot[i]);
    }
    printf(";");
    printf("\n=====================");
}



int genRandInt(int max)
/*
 * Purpose
 * =======
 * Generate a pseudo-random number in range [0, max];
 * It provides to avoid the bias problem;
 * Before invoking this function in a loop it's recommended to call srand();
 *
 * Arguments
 * =========
 * max: (int), upper limit of the range;
 *
 * Return
 * ======
 * (int) number in a range [0, max];
 */
{
    unsigned int
    numBins = (unsigned int) max + 1,
    numRand = (unsigned int) RAND_MAX + 1,
    binSize = numRand / numBins,
    defect   = numRand % numBins;

    int x;
    do {
        x = random();
    }
    while (numRand - defect <= (unsigned int)x);

    return x/binSize;
}



void freeMat(double **mat, int nRows)
/*
 * Purpose
 * =======
 * Deallocate memory of a matrix(double**);
 *
 * Arguments
 * =========
 * mat: (double**)[nRows,nRows], matrix to deallocate;
 * nRows: (int), number of rows of mat;
 */
{
    int i;
    for (i = 0; i < nRows; i++){
        free(mat[i]);
    }
    free(mat);
}


void freeMatInt(int **mat, int nRows)
/*
 * Purpose
 * =======
 * Deallocate memory of a matrix(int**);
 *
 * Arguments
 * =========
 * mat: (int**)[nRows,nRows], matrix to deallocate;
 * nRows: (int), number of rows of mat;
 */
{
    int i;
    for (i = 0; i < nRows; i++){
        free(mat[i]);
    }
    free(mat);
}


double** buildSample(char* filename, int dim)
/*
 * Purpose
 * =======
 * Build up a bootstrap replicate of a given dimension;
 * Currently not used.
 *
 * Arguments
 * =========
 * filename: (string), name of the dataset file;
 * dim: (int), size of the bootstrap replicate. It's equivalent at the number of rows to select up;
 *
 * Return
 * ======
 * matrix: (double**)[dim, colData], bootstrap replicate. This have the same format of the dataset except the number of rows;
 */
{
    int i, *indices;
    double **matrix;
    /* Build up an array with the indices of the rows to take */
    indices = malloc(sizeof(int)*dim);
    srand(time(NULL));
    for(i = 0; i<dim; i++){
        indices[i] = genRandInt(rowData);
    }
    /* Sort index */
    quicksort(indices, 0, dim-1);
    /* Get the matrix of the subset of the dataset*/
    matrix =  readSubsetCSVFile(filename, dim, colData, indices);
    free(indices);
    return matrix;
}



double** getMinor(double **distMat, int dim, int* indices)
/*
 * Purpose
 * =======
 * Get a minor of a matrix;
 * Note: currently not used.
 *
 * Arguments
 * =========
 * filename: (string), name of the dataset file;
 * dim: (int), size of the bootstrap replicate. It's equivalent at the number of rows to select up;
 * indices: (int*), array of size dim. Each element will be an index;
 *
 * Return
 * ======
 * matrix: (double**)[dim, dim], bootstrap replicate. This have the same format of the dataset except the number of rows;
 */
{
    int i, j;
    double **minorMat;
    /* Build up an array with the indices of the rows to take */
    srand(time(NULL));
    for(i = 0; i<dim; i++){
        indices[i] = genRandInt(rowData-1);
    }
    /* Sort index */
    quicksort(indices, 0, dim-1);

    if (dim < 2)
        return NULL;

    /* Set up the minor matrix */
    minorMat = malloc(dim*sizeof(double*));
    if(minorMat==NULL) return NULL; /* Not enough memory available */
    minorMat[0] = NULL; /* The zeroth row has zero columns. We allocate it anyway for convenience.*/

    for (i = 1; i < dim; i++)
    {
        minorMat[i] = malloc(i*sizeof(double));
        if (minorMat[i]==NULL)
            break; /* Not enough memory available */
    }
    if (i < dim) /* break condition encountered */
    {
        j = i;
        for (i = 1; i < j; i++)
            free(minorMat[i]);
        return NULL;
    }

    /* Calculate the distances and save them in the ragged array */
    for (i = 0; i < dim; i++){
        for (j = i+1; j < dim; j++){
            if (indices[i]==indices[j]){
                minorMat[j][i] = 0;
            }
            else{
                minorMat[j][i] = distMat[indices[j]][indices[i]];
            }
        }
     }
    return minorMat;
}



int** hc2split(Node* node, int nnodes, int toPrint)
/*
 * Purpose
 * =======
 * Extract clusters contents from the output of hiearchical clustering;
 *
 * Arguments
 * =========
 * node: (node), output of hierachical clustering;
 * nnodes: (int), number of nodes in the output of hierachical clustering;
 * toPrint: (int), if 1 enable the print of objects. Useful for testing;
 *
 * Return
 * ======
 * CC: (int**)[nnodes, rowData]. Matrix of 0 and 1. If C[i][j] is 1, the element j is in the cluster i;
 */
{
    int **CC;   //each row is a cluster and each column an element
    int B[nnodes][rowData]; //each row is a cluster. each column an element that belogs to the cluster
    int mergeM[nnodes][2];
    int i,j, ai, z, **A;

    /* Inizialize CC */
    CC = malloc(nnodes*sizeof(int*));
    for (i = 0; i<nnodes; i++){
        CC[i] = malloc((rowData)*sizeof(int*));
    }
    for (i = 0; i < nnodes; i++)
        for (j = 0; j < rowData; j++)
            CC[i][j] = 0;

    /* Fill MergeM */
    /* MergeM describes the output of clustering. It is matrix[nnodes, 2].
     * Each row contains tha distance between two objects.
     * If the element is positive it's a single element, else it's a cluster.
     * */
    for(i=0; i<nnodes; i++){
        mergeM[i][0] = node[i].left;
        mergeM[i][1] = node[i].right;
    }

    /* Inizialize B */
    for (i = 0; i < nnodes; i++)
        for (j = 0; j < rowData; j++)
            B[i][j] = -1;

    /*
     * Fill B in the proper way.
     * Rows of B are cluster, while columns are elements.
     */
    for(i = 0; i < nnodes; i++){
        ai = mergeM[i][0];
        j=0;
        if (ai >= 0){
            B[i][j] = ai;
        }
        else if (ai < 0){
            ai = -ai;
            while(B[ai-1][j] >= 0){
                B[i][j] = B[ai-1][j];
                j++;
            }
        }
        ai = mergeM[i][1];
        if (j==0) j++;
        if (ai >= 0)
            B[i][j] = ai;
        else if (ai < 0){
            ai = -ai;
            z=j;
            j=0;
            while(B[ai-1][j] >= 0){
                B[i][z] = B[ai-1][j];
                j++;
                z++;
            }
        }
    } //end for i

    /* Replace -1 values with MAX_INT for sorting */
    for (i=0; i<nnodes;i++)
        for (j=0; j<rowData; j++)
            if (B[i][j]<0)
                B[i][j] = INT_MAX;
    /* Sort row of B */
    for (i=0; i<nnodes;i++)
        quicksort(B[i],0,nnodes-1);
    /* Fill CC*/
    for(i = 0; i < nnodes; i++){
        for(j = 0; j < rowData; j++)
            if(B[i][j] != INT_MAX){
                CC[i][B[i][j]] = 1;
            }
    }

    /* Print */
    if (toPrint == 1){
        A = malloc(nnodes*sizeof(int*));
        for (i = 0; i<nnodes; i++){
            A[i] = malloc((rowData)*sizeof(int*));
        }
        for (i = 0; i < nnodes; i++)
            for (j = 0; j < rowData; j++)
                A[i][j] = B[i][j];
        printf("\n==== A ====");
        showDataInt(nnodes, rowData, A);
        printf("==== CC ====");
        showDataInt(nnodes, rowData, CC);
    }

    return CC;
}



int** buildPattern()
/*
 * Purpose
 * =======
 * Get the output of clustering on the complete data.
 *
 * It's been called from R for make comparisons.
 *
 * Arguments
 * =========
 *
 * Return
 * ======
 * CC: (int**)[rowData-1, rowData]. Matrix of 0 and 1. If C[i][j] is 1, the element j is in the cluster i;
 */
{
    int j, k, **globalC, **mask=malloc(rowData*sizeof(int*));
    double  **dataMat, **distance, **mat;
    char typeDist = 'c', method = 'a';
    Node* tree;

    /* Calculate global distance matrix */
    dataMat = readCSVFile(filen);
    /* Tranpose */
    mat = transposeMat(dataMat, colData, rowData);

    /* Hierachical clustering of complete data */
    distance = compDistMatrix(rowData, colData, mat, 0);
    for (k=0;k<rowData;k++)
        mask[k] = malloc(colData*sizeof(int));
    for (k=0;k<rowData;k++)
        for(j=0;j<colData;j++)
            mask[k][j]=1;
    tree = hclusteringM(rowData, colData, mat, mask, distance, typeDist, method, 0);
    freeMat(dataMat, colData);
    globalC=hc2split(tree, rowData-1, 0);

    return globalC;
}



void updateFrequency(int** globalC, int** localC, int nClust,int nRows, int *freq){
/*
 * Purpose
 * =======
 * Update frequency values of freq by matching elements in common between globalC and localC;
 * It writes in freq;
 *
 * Arguments
 * =========
 * globalC: (int**)[nClust,rowData], set of clusters of the global tree;
 * localC: (int**)[nClust,rowData], set of clusters of the local tree (bootstrap replicate);
 * nClust: (int), number of clusters;
 * nRows: (int), number of elements;
 * freq: (int*)[nClust], array with frequencies to update;
 *
 */
    int i,j, z, flag;

    for (i=0;i<nClust; i++){    //for each row of localB
        for (j=0;j<nClust;j++){ //for each row of globalC
            flag=0;
            z=0;
            while (flag==0 && z<nRows){  //compare each element of the two rows
                if (localC[i][z] != globalC[j][z])
                    flag=1;
                z++;
            }
            if (flag==0)
                freq[j]++;
        }   //end for j
     }  //end for i
}


int** eboot(int toPrint)
/*
 * Pvclust main function;
 *
 * First of all it does an hierarchical clustering on complete data (all objects, all features).
 * Then it does multiscale bootstrap. For each values in lenDimsBoot computes nBoot iterations.
 * In each of those iteration is performed a resampling with hierachical clustering.
 * It's calculated how many times a cluster belonging at tree of the complete data appears in multiscale bootstrap.
 *
 * Arguments
 * =========
 * toPrint: (int), if 1 print the frequencies
 *
 * Return
 * ======
 * freq: (int**)[lenDimsBoot, rowData-1], freq[i][j] says the number of times the cluster j has been occurred in the replication i.
 */
{
    int i,j, k, nCols, *indices, **globalC, **localC, **freq, **mask=malloc(rowData*sizeof(int*));
    double **bootSample, **dataMat, **distance, **mat, **BP;
    char typeDist = 'd', method = 'a';
    Node* tree;

    /* Calculate global distance matrix */
    dataMat = readCSVFile(filen);
    /* Tranpose */
    mat = transposeMat(dataMat, colData, rowData);

    /* Hierachical clustering of complete data */
    distance = compDistMatrix(rowData, colData, mat, 0);
    for (k=0;k<rowData;k++)
        mask[k] = malloc(colData*sizeof(int));
    for (k=0;k<rowData;k++)
        for(j=0;j<colData;j++)
            mask[k][j]=1;

    /* Global tree */
    tree = hclusteringM(rowData, colData, mat, mask, distance, typeDist, method, 0);
    freeMat(dataMat, colData);
    globalC=hc2split(tree, rowData-1, 0);

    /* freq is a matrix of 0 at the start.
     * freq is associate at globalC:
     * freq[i][j] stands for the frequency of j-th cluster described by j-th row of globalCC
     * i stands for the i-th size of the bootstrap
     */
    freq = malloc((lenDimsBoot)*sizeof(int*));
    for (i=0;i<lenDimsBoot;i++){
        freq[i] = malloc((rowData-1)*sizeof(int));
    }
    for (i=0;i<lenDimsBoot;i++)
        for (j=0;j<rowData-1;j++)
            freq[i][j]=0;

    /* Bootstrap */
    for(i=0; i<lenDimsBoot; i++){   //for each size of bootstrap
        nCols = (int)(colData * dimsBoot[i])/100;   //it's the size of the sample of the boot (number of columns)
        printf("\n--> Bootstrap (r = %.2f, nCols = %d);",((double)dimsBoot[i]/100), nCols);
        indices = malloc(sizeof(int)*nCols);
        /* missing data mask for hclusteringM*/
        for (k=0;k<rowData;k++)
            mask[k] = malloc(nCols*sizeof(int));
        for (k=0;k<rowData;k++)
            for(j=0;j<nCols;j++)
                mask[i][j]=1;
        /* for each bootstrap replications */
        for (j=0; j<nBoot; j++){
            //printf("\nj: %d", j);
            /* Build up the sample */
            for(k = 0; k<nCols; k++){
                srand ( time(NULL)+rand());
                indices[k] = genRandInt(colData-1); //sampling
            }
            quicksort(indices, 0, nCols-1);
            bootSample = loadSubsetDataCol(mat, rowData, nCols, indices);
            /* Compute distance matrix */
            distance = compDistMatrix(rowData, nCols, bootSample, 0);
            /* Hierarchical clustering */
            tree = hclusteringM(rowData, nCols, bootSample, mask, distance, typeDist, method, 0);
            /* Get cluster frequency */
            localC=hc2split(tree, rowData-1, 0);
            /* Update frequency */
            updateFrequency(globalC, localC, rowData-1, rowData, freq[i]);
            /* free */
            freeMat(bootSample, rowData);   // or memory leak :)
            freeMatInt(localC, rowData-1);
            free(tree);
            freeMat(distance, rowData);
        } //end for j
        free(indices);
        for(k = 0; k<rowData; k++)  free(mask[k]);
    } //end for i
    freeMat(mat, rowData);

    /* Print */
    if (toPrint == 1){
        printf("\n\n===== Frequency =====");
        for (i=0;i<lenDimsBoot;i++){
            printf("\nR->");
            for (j=0;j<rowData-1;j++){
                printf("\t%d",freq[i][j]);
            }
        }
    }//end print

    return freq;
}



void epvclust(char* filename, int nboot)
/*
 * Purpose
 * =======
 * Run pvclust
 *
 * Arguments
 * =========
 * filname: (char*), name of file with data. It must have not names on columns and rows and need to be in CSV format
 *          (and in particular data must be numerics and separated by comma)
 * nboot: (int), number of times of bootstrap interation for each value of r
 */
{
    int** m;

    /* Inizialize variables */
    pvclustInit(filename, nboot, 0, 0);
    /* Pvclust */
    m = eboot(0);

    return;
}


int main(void)
/*
 * The program works with datasets that have objects on the columns and feature on the rows.
 * i.e. patiens on the columns and genes on the rows.
 *
 * Later in eboot() the matrix will be transpose for convenience of convention.
 *
 * This is an example of execution.
 */
{
    int nboot = 100;
    char *filename = "./../exampleData/medium.csv"; //path dataset
    clock_t begin, end;
    double time_spent;

    printf("Starting..");

    begin = clock();

    printf("\n");
    epvclust(filename, nboot);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("\n\nTime elapsed: %f", time_spent);

    printf("\nEnding..\n");

    return 0;
}
