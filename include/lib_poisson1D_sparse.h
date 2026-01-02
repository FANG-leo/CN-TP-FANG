#ifndef LIB_POISSON1D_SPARSE_H
#define LIB_POISSON1D_SPARSE_H

void poisson1D_CSR(int n, double **val, int **col, int **rowptr);
void dcsrmv(int n, double *val, int *col, int *rowptr,
            double *x, double *y);

void poisson1D_CSC(int n, double **val, int **row, int **colptr);
void dcscmv(int n, double *val, int *row, int *colptr,
            double *x, double *y);

#endif