/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator
  int n = *la;
  int lab_v = *lab;
  int ku = *kv; if (ku <= 0) ku = 1;
  double h = 1.0/(n+1.0);
  double inv_h2 = 1.0/(h*h);
  int size = lab_v * n;
  for (int k = 0; k < size; ++k) AB[k] = 0.0;
  for (int j = 0; j < n; ++j) {
    int row_diag = ku + j - j; /* = ku */
    AB[row_diag + j*lab_v] = 2.0 * inv_h2;
    if (j+1 < n) AB[(ku + 1) + j*lab_v] = -1.0 * inv_h2;
    if (j-1 >= 0) AB[(ku - 1) + j*lab_v] = -1.0 * inv_h2;
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
  int n = *la;
  int lab_v = *lab;
  int ku = *kv;
  if (ku <= 0) ku = 1;

  int size = lab_v * n;
  for (int i = 0; i < size; ++i) AB[i] = 0.0;

  for (int j = 0; j < n; ++j) {
    AB[ku + j*lab_v] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
  int n = *la;
  double h = 1.0 / (n + 1.0);
  double inv_h2 = 1.0 / (h*h);

  for (int i = 0; i < n; ++i) RHS[i] = 0.0;
  RHS[0] += (*BC0) * inv_h2;
  RHS[n-1] += (*BC1) * inv_h2;
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
  int n = *la;
  double T0 = *BC0, T1 = *BC1;
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
    num += d*d;
    den += y[i]*y[i];
  }
  if (den == 0.0) return (num==0.0) ? 0.0 : 1e300;
  return sqrt(num)/sqrt(den);
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  int ku = 1; /* tri-diagonal case */
  int lab_v = *lab;
  int row = ku + i - j;
  return row + j * lab_v;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  int nn = *la;           /* actual matrix size */
  (void)n;                /* keep signature but mark unused to avoid warnings */

  int lab_v = *lab;
  int ku_v = *ku;
  if (ku_v <= 0) ku_v = 1;

  /* initialize ipiv to identity (no pivoting) */
  for (int i = 0; i < nn; ++i) ipiv[i] = i+1;

  *info = 0;
  double eps = 1e-18;

  /* Thomas-like elimination stored in GB layout */
  for (int j = 0; j < nn-1; ++j) {
    /* pivot at A(j,j) */
    int idx_pivot = ku_v + j * lab_v; /* row = ku, col = j */
    double pivot = AB[idx_pivot];

    if (fabs(pivot) <= eps) {
      *info = j+1; /* Fortran-style positive index of failure */
      return *info;
    }

    /* sub-diagonal A(j+1,j) stored at row = ku+1 in column j */
    int idx_sub = (ku_v + 1) + j * lab_v;
    double mult = AB[idx_sub] / pivot;
    AB[idx_sub] = mult; /* store multiplier (L) */

    /* upper element A(j,j+1) is in column j+1 at row = ku - 1 */
    int idx_upper_in_nextcol = (ku_v - 1) + (j+1) * lab_v;
    /* diagonal at next column */
    int idx_diag_next = ku_v + (j+1) * lab_v;

    /* update next diagonal */
    AB[idx_diag_next] -= mult * AB[idx_upper_in_nextcol];
  }

  *info = 0;
  return *info;
}
