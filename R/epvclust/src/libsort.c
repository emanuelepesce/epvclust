/* 
* Sorting algorithms library
*/
#include <stdio.h>
#include <stdlib.h>
#include "libsort.h"

//=============================================================
/*
* Quicksort algorithm. Code by "Rosetta code"
*/
void swap_r(int r[], int a, int b)
{
    int temp = r[a];
    r[a] = r[b];
    r[b] = temp;
}
 
void quicksort(int r[], int start, int end)
{
    if(end > start)
    {
        int pivot_index = (start + end) / 2;
        int pivot = r[pivot_index];
        int chg, i;
 
        swap_r(r, pivot_index, end);
 
        for(i = chg = start; i < end; i++)
        {
            if(r[i] < pivot)
            {
                swap_r(r, i, chg);
                chg++;
            }
        }
 
        swap_r(r, chg, end);
 
        quicksort(r, start, chg - 1);
        quicksort(r, chg + 1, end);
    }
}

//=============================================================
