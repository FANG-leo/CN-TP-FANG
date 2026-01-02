#include <stdlib.h>
#include "lib_poisson1D_sparse.h"


void poisson1D_CSR(int n, double **val, int **col, int **rowptr)
{
  int nnz = 3*n - 2;

  *val    = (double*) malloc(nnz * sizeof(double));
  *col    = (int*)    malloc(nnz * sizeof(int));
  *rowptr = (int*)    malloc((n+1) * sizeof(int));

  int k = 0;
  (*rowptr)[0] = 0;

  for (int i = 0; i < n; i++) {

    if (i > 0) {
      (*val)[k] = -1.0;
      (*col)[k] = i - 1;
      k++;
    }

    (*val)[k] = 2.0;
    (*col)[k] = i;
    k++;

    if (i < n-1) {
      (*val)[k] = -1.0;
      (*col)[k] = i + 1;
      k++;
    }

    (*rowptr)[i+1] = k;
  }
}

void dcsrmv(int n, double *val, int *col, int *rowptr,
            double *x, double *y)
{
  for (int i = 0; i < n; i++) {
    y[i] = 0.0;
    for (int k = rowptr[i]; k < rowptr[i+1]; k++) {
      y[i] += val[k] * x[col[k]];
    }
  }
}


void poisson1D_CSC(int n, double **val, int **row, int **colptr)
{
  int nnz = 3*n - 2;

  *val    = (double*) malloc(nnz * sizeof(double));
  *row    = (int*)    malloc(nnz * sizeof(int));
  *colptr = (int*)    malloc((n+1) * sizeof(int));

  int k = 0;
  (*colptr)[0] = 0;

  for (int j = 0; j < n; j++) {

    if (j > 0) {
      (*val)[k] = -1.0;
      (*row)[k] = j - 1;
      k++;
    }

    (*val)[k] = 2.0;
    (*row)[k] = j;
    k++;

    if (j < n-1) {
      (*val)[k] = -1.0;
      (*row)[k] = j + 1;
      k++;
    }

    (*colptr)[j+1] = k;
  }
}


void dcscmv(int n, double *val, int *row, int *colptr,
            double *x, double *y)
{
  for (int i = 0; i < n; i++)
    y[i] = 0.0;

  for (int j = 0; j < n; j++) {
    for (int k = colptr[j]; k < colptr[j+1]; k++) {
      y[row[k]] += val[k] * x[j];
    }
  }
}
