/*
* libwwd.h
*/
#ifndef LIBWWD_H_
#define LIBWWD_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

int getNumberOfCols(char* line, char* del);
int getNumberOfRows(FILE *file);
char *str_replace (const char *string, const char *substr, const char *replacement );
double** loadData(FILE *file, int nCols, int nRows, char *del);
double** readCSVFile(char* filename);
int* getDatasetSize(char *filename);
void showData(int nRows, int nCols, double** data);
double** loadSubsetData(FILE *file, int nCols, char *del, int *indices, int sizeI);
double** readSubsetCSVFile(char* filename, int nRows, int nCols, int *indices);
void showLowerMatrix(int nRows, double **data);
void showDataInt(int nRows, int nCols, int** data);
double** readSubsetCSVFileCol(char* filename, int nRows, int nCols, int *indices);
double** loadSubsetDataCol(double** cmatrix, int nRows, int nCols, int *indices);
double** transposeMat(double** mat, int nRows, int nCols);


#endif /*libwwd.h*/


