/*
 *	A collection of useful functions for working with datasets
 *  
 *  Author: Emanuele Pesce
 *
 *  Last Update: 10 December 2015
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "libwwd.h"


double** readCSVFile(char* filename){
/*
    Purpose
    =======
    Load in a matrix of double a CSV dataset.
	It works if csv file has only numerics values.
	
    arguments
    =========
    filename: char*, string that cointains name of the dataset file

    return
    ======
	matrix: double**, matrix of double that cointains dataset values
*/
	double** matrix;
	char line[BUFSIZ];
	int nCols, nRows;
	FILE *fp;
	int *dim;
	
	/*	Open the file	*/
	fp = fopen(filename, "r");
	if (fp != NULL){
		/*	Get the number of columns and rows of the dataset	*/
		dim=getDatasetSize(filename);
		nRows=dim[0];
		nCols=dim[1];

		/*	Load dataset in a matrix of double	*/
		matrix = loadData(fp, nRows, nCols, ",");
	}
	/*error with the file opening */
	else{ 
		perror(filename);
	}
	return matrix;
}



double** loadData(FILE *file, int nRows, int nCols, char *del){
/*	
    Purpose
    =======
	Put the dataset values in a matrix.
	It assumes that dataset cointains only numerics values.	

    arguments
    =========
    file: FILE, file pointer of the dataset
    nCols: int, number of columns of the dataset
    nRows: number of rows of the dataset

    return
    ======
	matrix: double**, a matrix that cointains the dataset values
*/
	int i,j; 
	char line[BUFSIZ], *stmp;
	double** matrix = malloc(nRows*sizeof(double*) );
	
	/*	Allocate matrix */
	for (i = 0; i < nRows; i++){
		matrix[i] = malloc(nCols*sizeof(double));
    }
	/*
	  	For each line, remove the separators and put elements in a matrix as
	 	double
	 */
	rewind(file);
	for ( i = 0; fgets(line, sizeof line, file); ++i ){
	     stmp = line;
	     stmp=str_replace(stmp, del, " ");
	     for ( j = 0; j<nCols; j++)
	     {		     
	        matrix[i][j] = strtold(stmp, &stmp);
	     }
	  }
	fclose(file);
	return matrix;
}



int getNumberOfCols(char* line, char* del){
/*
    Purpose
    =======
	Count the number of columns of a row.
	
    arguments
    =========
    line: string, row of a dataset
    del: string, separator of the string values of line

    return
    ======
	nCols: number of columns of a row
*/
	char* col;
	int nCols = 0;
	
	for (col = strtok(line,del); col != NULL; col = strtok(NULL, del)){
		nCols++;
	}
	
	return nCols;
}


int getNumberOfRows(FILE *file){
/*
    Purpose
    =======
    Count the number of rows of a file
	It assumes that the last line does not cointan "\n"
	
	arguments
    =========
	file: FILE, file pointer

    return
    ======
    nRows: number of rows of a file
*/
	int cur, nRows=0;
	
	while ( (cur=fgetc(file)) != EOF ) {
        if ( cur == '\n' )
            nRows++;
    }
	return nRows;
}



char *str_replace (const char *string, const char *substr, const char *replacement ){
/*
    Purpose
    =======
	implements a str_replace PHP like function	
	
	arguments
    =========
	string: string, string to modify
	substring: string, substring to change in string
	replacement: string, string to substitute in string at place of substring

    return
    ======
    newstring: string, string modified
*/
  char *tok = NULL;
  char *newstr = NULL;
  char *oldstr = NULL;
  char *head = NULL;
 
  /* if either substr or replacement is NULL, duplicate string a let caller handle it */
  if ( substr == NULL || replacement == NULL ) return strdup (string);
  newstr = strdup (string);
  head = newstr;
  while ( (tok = strstr ( head, substr ))){
    oldstr = newstr;
    newstr = malloc ( strlen ( oldstr ) - strlen ( substr ) + strlen ( replacement ) + 1 );
    /*failed to alloc mem, free old string and return NULL */
    if ( newstr == NULL ){
      free (oldstr);
      return NULL;
    }
    memcpy ( newstr, oldstr, tok - oldstr );
    memcpy ( newstr + (tok - oldstr), replacement, strlen ( replacement ) );
    memcpy ( newstr + (tok - oldstr) + strlen( replacement ), tok + strlen ( substr ), strlen ( oldstr ) - strlen ( substr ) - ( tok - oldstr ) );
    memset ( newstr + strlen ( oldstr ) - strlen ( substr ) + strlen ( replacement ) , 0, 1 );
    /* move back head right after the last replacement */
    head = newstr + (tok - oldstr) + strlen( replacement );
    free (oldstr);
  }
  return newstr;
}


int* getDatasetSize(char *filename){
/*
    Purpose
    =======
	Counts the number of columns and rows of a dataset file
	
	arguments
    =========
	filename: char*, string containing the name of the file

    return
    ======
    dim: int*, dim[1]:number of rows, dim[2]: number of columns
*/
	int *dim;
	char line[BUFSIZ];
	FILE *fp;
	
	dim = malloc(2*sizeof(int));
	
	/* Open the file	*/
	fp = fopen(filename, "r");
	if (fp != NULL){
		/*	Get the number of columns and rows of the dataset	*/
		dim[0]=getNumberOfRows(fp);
		rewind(fp);
		fgets(line, sizeof line, fp);
		dim[1]=getNumberOfCols(line, ",");
		
	}
	/*error with the file opening */
	else{ 
		perror(filename);
	}
	fclose(fp);
	return dim;
}


void showData(int nRows, int nCols, double** data){
/*
    Purpose
    =======
	Print the data matrix
	
	arguments
    =========
	nrows: int number of rows
	ncols: int, number of columns
	data: double **, matrix of data

*/
	int i, j;	
	printf("\n=========== Data =============\n");
  	for (j = 0; j < nCols; j++) 
  		printf("\tCol %d", j);
  	printf ("\n");
  	for (i = 0; i < nRows; i++){
  	   	printf("Row %d:", i);
    	for (j = 0; j < nCols; j++){ 
    	  	printf("\t%f",data[i][j]);
    	}
    printf("\n");
  }
  printf("\n");
  return;
}

void showDataInt(int nRows, int nCols, int** data){
/*
    Purpose
    =======
    Print the data matrix

    arguments
    =========
    nrows: int number of rows
    ncols: int, number of columns
    data: double **, matrix of data
*/
    int i, j;
    printf("\n=========== Data =============\n");

    for (j = 0; j < nCols; j++)
        printf("\tCol %d", j);
    printf ("\n");
    for (i = 0; i < nRows; i++){
        printf("Row %d:", i);
        for (j = 0; j < nCols; j++){
            printf("\t%d",data[i][j]);
        }
    printf("\n");
  }
  printf("\n");
  return;
}


double** loadSubsetData(FILE *file, int nCols, char *del, int *indices, int sizeI){
/*
    Purpose
    =======
    Put a subset dataset values in a matrix.
    It takes care only about the rows specified in indices.
    It assumes that dataset cointains only numerics values.

    arguments
    =========
    file: FILE, file pointer of the dataset
    nCols: int, number of columns of the dataset
    del: char*, delimitator
    indices: int*, array containing indices to select. It MUST BE sorted.
    sizeI: int, lenght of indices (alias number of rows of the subset dataset)

    return
    ======
    matrix: double**, a matrix that cointains the subset dataset values
*/
    int i,j,numLine = 0;
    char line[BUFSIZ], *stmp;
    double** matrix = malloc(sizeI*sizeof(double*) );

    /*	Allocate matrix */
    for (i = 0; i < sizeI; i++){
        matrix[i] = malloc(nCols*sizeof(double));
    }
    /*
        For each line to select, remove the separators and put elements in a matrix as
        double
     */
    rewind(file);
    fgets(line, sizeof line, file);
    i=0;
    while (i < sizeI){
        /* If the index is the same increase "i" else go on with next line of the file */
        if (indices[i] == numLine){
            stmp = line;
            stmp=str_replace(stmp, del, " ");
            for ( j = 0; j<nCols; j++)
            {
               matrix[i][j] = strtold(stmp, &stmp);
            }
            i++;
        }
        else{
            fgets(line, sizeof line, file);
            numLine++;
        }
    }
    fclose(file);
    return matrix;
}

double** readSubsetCSVFileCol(char* filename, int nRows, int nCols, int *indices){
/*
    Purpose
    =======
    Load in a matrix of double a a part CSV dataset.
    It works if csv file has only numerics values.
    All the rows will be taken, but not all the columns.


    arguments
    =========
    filename: char*, string that cointains name of the dataset file
    nRows: int, number of rows.
    nCols: int, number of columns. It's equal to the number of elements of indices.
    indices: int*, array. Each element is an index a row to select up.

    return
    ======
    matrix: double**[nRows,nCols], matrix of double that cointains dataset values.

*/
    double **cmatrix, **matrix;
    int i, j=0;
    FILE *fp;
    /*	Open the file	*/
    /*	Allocate matrix */
    matrix = malloc(nRows*sizeof(double*)) ;
    for (i = 0; i < nRows; i++){
        matrix[i] = malloc(nCols*sizeof(double));
    }
    fp = fopen(filename, "r");
    if (fp != NULL){
        /*	Load dataset in a matrix of double	*/
        cmatrix = readCSVFile(filename);

        /* Fill matrix */
        for (j=0;j<nCols;j++){
            for (i=0;i<nRows;i++){
                matrix[i][j] = cmatrix[i][indices[j]];
            }
        }
    }
    /*error with the file opening */
    else{
        perror(filename);
    }
    return matrix;
}



double** readSubsetCSVFile(char* filename, int nRows, int nCols, int *indices){
/*
    Purpose
    =======
    Load in a matrix of double a CSV dataset.
    It works if csv file has only numerics values.

    arguments
    =========
    filename: char*, string that cointains name of the dataset file
    nRows: int, number of rows. It's equal to the number of elements of indices.
    nCols: int, number of columns
    indices: int*, array. Each element is an index a row to select up.

    return
    ======
    matrix: double**, matrix of double that cointains dataset values
*/
    double** matrix;
    FILE *fp;
    /*	Open the file	*/
    fp = fopen(filename, "r");
    if (fp != NULL){

        /*	Load dataset in a matrix of double	*/
        matrix = loadSubsetData(fp, nCols, ",", indices, nRows);
    }
    /*error with the file opening */
    else{
        perror(filename);
    }
    return matrix;
}


double** loadSubsetDataCol(double** cmatrix, int nRows, int nCols, int *indices){
/*
    Purpose
    =======
    Load in a matrix, a subset of cmatrix.
    All rows will be taken. It's not the same for columns because only those specified in
    indices will be taken.

    arguments
    =========
    filename: cmatrix, complete matrix of data.
    nRows: int, number of rows.
    nCols: int, number of columns. It's equal to the number of elements of indices.
    indices: int*, array. Each element is an index a row to select up.

    return
    ======
    matrix: double**[nRows,nCols], matrix of double that cointains dataset values.
*/
    double** matrix;
    int i, j=0;
    /*	Allocate matrix */
    matrix = malloc(nRows*sizeof(double*) );
    for (i = 0; i < nRows; i++){
        matrix[i] = malloc(nCols*sizeof(double));
    }
    /* Fill matrix */
    for (j=0;j<nCols;j++){
        for (i=0;i<nRows;i++){
            matrix[i][j] = cmatrix[i][indices[j]];
        }
    }

    return matrix;
}


void showLowerMatrix(int nRows, double **data)
/*
    Purpose
    =======
    Print the matrix "data". It assumes that data is a lower triangular matrix;

    arguments
    =========
    data: double**, lower traingular matrix of double;
    nRows: int, number of rows of data;
*/
{
    int i,j;

    printf("\n============ Matrix ============\n");
    printf("   Obj:");
    for(i=0; i<nRows-1; i++)
        printf("%6d", i);
    printf("\n");
    for(i=0; i<nRows; i++)
    {
        printf("Obj %2d:",i);
        for(j=0; j<i; j++) printf(" %5.2f",data[i][j]);
        printf("\n");
    }
    printf("\n");

}


double** transposeMat(double** mat, int nRows, int nCols){
/*
    Purpose
    =======
    Transpose a matrix

    return
    ======
    tmat, double**, transpose of mat
*/
    double** tmat;
    int i,j;

    tmat = malloc(nCols*sizeof(double*));
    for (i = 0; i < nCols; i++){
        tmat[i] = malloc(nRows*sizeof(double));
    }
    for (i=0;i<nRows;i++){
        for (j=0;j<nCols;j++){
            tmat[j][i]=mat[i][j];
        }
    }
    return tmat;
}
