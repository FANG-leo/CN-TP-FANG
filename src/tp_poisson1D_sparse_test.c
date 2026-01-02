#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#include "lib_poisson1D.h"
#include "lib_poisson1D_sparse.h"

int main(void)
{
  int n = 10;

  /* ---------- Dense vector x ---------- */
  double *x = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++)
    x[i] = 1.0;   /* simple test vector */

  /* ---------- GB matrix ---------- */
  int kv = 0, ku = 1, kl = 1;
  int lab = kv + ku + kl + 1;

  double *AB = (double*) malloc(lab * n * sizeof(double));
  set_GB_operator_colMajor_poisson1D(AB, &lab, &n, &kv);

  double *y_gb = (double*) calloc(n, sizeof(double));

  cblas_dgbmv(CblasColMajor, CblasNoTrans,
              n, n, kl, ku,
              1.0, AB, lab,
              x, 1, 0.0, y_gb, 1);

  /* ---------- CSR matrix ---------- */
  double *val_csr;
  int *col_csr, *rowptr;

  poisson1D_CSR(n, &val_csr, &col_csr, &rowptr);

  double *y_csr = (double*) malloc(n * sizeof(double));
  dcsrmv(n, val_csr, col_csr, rowptr, x, y_csr);

  /* ---------- CSC matrix ---------- */
  double *val_csc;
  int *row_csc, *colptr;

  poisson1D_CSC(n, &val_csc, &row_csc, &colptr);

  double *y_csc = (double*) malloc(n * sizeof(double));
  dcscmv(n, val_csc, row_csc, colptr, x, y_csc);

  /* ---------- Compare results ---------- */
  double err_csr = 0.0, err_csc = 0.0;
  for (int i = 0; i < n; i++) {
    err_csr += (y_gb[i] - y_csr[i]) * (y_gb[i] - y_csr[i]);
    err_csc += (y_gb[i] - y_csc[i]) * (y_gb[i] - y_csc[i]);
  }

  err_csr = sqrt(err_csr);
  err_csc = sqrt(err_csc);

  printf("||Ax_GB - Ax_CSR||_2 = %e\n", err_csr);
  printf("||Ax_GB - Ax_CSC||_2 = %e\n", err_csc);

  /* ---------- Cleanup ---------- */
  free(x);
  free(AB);
  free(y_gb);
  free(val_csr);
  free(col_csr);
  free(rowptr);
  free(y_csr);
  free(val_csc);
  free(row_csc);
  free(colptr);
  free(y_csc);

  return 0;
}
