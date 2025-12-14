/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator
  int n    = *la;         
  int ldab = *lab;         

  double h    = 1.0 / (n + 1);
  double diag =  2.0 / (h * h);
  double off  = -1.0 / (h * h);

  for (int j = 0; j < n; j++) {
    for (int i = 0; i < ldab; i++) {
      AB[indexABCol(i, j, lab)] = 0.0;
    }
  }

  int row_extra  = 0; 
  (void)row_extra;   
  int row_super  = 1;  
  int row_diag   = 2;  
  int row_sub    = 3;  

  for (int j = 0; j < n; j++) {
    if (j > 0) {
      AB[indexABCol(row_super, j, lab)] = off;
    }

    AB[indexABCol(row_diag, j, lab)] = diag;

    if (j < n - 1) {
      AB[indexABCol(row_sub, j, lab)] = off;
    }
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
  int n = *la;
  int ldab = *lab;
  int kv_v = *kv;
  int ku = kv_v; if (ku <= 0) ku = 1;

  int size = ldab * n;
  for (int k = 0; k < size; ++k) AB[k] = 0.0;
  for (int j = 0; j < n; ++j) {
    AB[indexABCol(ku, j, lab)] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
  int n = *la;
  double h = 1.0 / (n + 1.0);
  double inv_h2 = 1.0 / (h*h);
  for (int i = 0; i < n; ++i) RHS[i] = 0.0;
  RHS[0]     += (*BC0) * inv_h2;
  RHS[n-1]   += (*BC1) * inv_h2;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
  int n = *la;
  double T0 = *BC0;
  double T1 = *BC1;
  for (int i = 0; i < n; ++i) EX_SOL[i] = T0 + X[i] * (T1 - T0);
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
  int n = *la;
  double h = 1.0 / (n + 1.0);
  for (int i = 0; i < n; ++i) x[i] = (i + 1) * h;
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  int n = *la;
  double num = 0.0, den = 0.0;
  for (int i = 0; i < n; ++i){
    double d = x[i] - y[i];
    num += d * d;
    den += y[i] * y[i];
  }
  if (den == 0.0) {
    if (num == 0.0) return 0.0;
    return 1e300;
  }
  return sqrt(num) / sqrt(den);
}

int indexABCol(int row, int col, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  int ldab = *lab;
  return row + col * ldab;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  int nn = *la;
  (void)n;

  int lab_v = *lab;
  int kl_v  = *kl;
  int row_diag = lab_v - kl_v - 1;
  if (row_diag < 0) row_diag = 0;


  for (int i = 0; i < nn; ++i) ipiv[i] = i + 1;

  *info = 0;
  double eps = 1e-18;

  for (int j = 0; j < nn - 1; ++j) {
    int idx_pivot = row_diag + j * lab_v;
    double pivot = AB[idx_pivot];
    if (fabs(pivot) <= eps) {
      *info = j + 1; 
      return *info;
    }

    int idx_sub = (row_diag + 1) + j * lab_v;
    double mult = AB[idx_sub] / pivot;
    AB[idx_sub] = mult;

    int idx_upper_in_nextcol = (row_diag - 1) + (j + 1) * lab_v;
    int idx_diag_next = row_diag + (j + 1) * lab_v;

    AB[idx_diag_next] -= mult * AB[idx_upper_in_nextcol];
  }

  *info = 0;
  return *info;
}
