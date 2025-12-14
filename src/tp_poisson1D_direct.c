/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/* using direct methods (LU factorization)*/
/******************************************/
#include "lib_poisson1D.h"
#include <time.h>

#define TRF 0  /* Use LAPACK dgbtrf for LU factorization */
#define TRI 1  /* Use custom tridiagonal LU factorization */
#define SV 2   /* Use LAPACK dgbsv (all-in-one solver) */

/**
 * Main function to solve the 1D Poisson equation -u''(x) = f(x) with Dirichlet BC.
 * 
 * @param argc: Number of command-line arguments
 * @param argv: Array of argument strings
 *              argv[1] (optional): Implementation method (0=TRF, 1=TRI, 2=SV)
 * @return 0 on success
 */
int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;                      /* Error code for various operations */
  int jj;                        /* Loop counter */
  int nbpoints, la;              /* nbpoints: total points, la: interior points */
  int ku, kl, kv, lab;           /* Band matrix parameters: ku/kl=super/sub diagonals, kv=extra space, lab=leading dimension */
  int *ipiv;                     /* Pivot indices for LU factorization */
  int info = 1;                  /* LAPACK info parameter (0=success) */
  int NRHS;                      /* Number of right-hand sides */
  int IMPLEM = 0;                /* Implementation method (TRF, TRI, or SV) */
  double T0, T1;                 /* Boundary conditions: T0 at x=0, T1 at x=1 */
  double *RHS, *EX_SOL, *X;      /* RHS: right-hand side, EX_SOL: exact solution, X: grid points */
  double **AAB;                  /* Unused variable */
  double *AB;                    /* Coefficient matrix in band storage */

  double relres;                 /* Relative forward error */

  struct timespec t_start, t_end;
  double elapsed;


  if (argc == 2) {
    IMPLEM = atoi(argv[1]);
  } else if (argc > 2) {
    perror("Application takes at most one argument");
    exit(1);
  }

  /* Problem setup */
  NRHS=1;           /* Solving Ax=b with one right-hand side */
  nbpoints=1000;      /* Total number of discretization points (including boundaries) */
  la=nbpoints-2;    /* Number of interior points (excluding boundaries) */
  T0=-5.0;          /* Dirichlet boundary condition at x=0 */
  T1=5.0;           /* Dirichlet boundary condition at x=1 */

  printf("--------- Poisson 1D ---------\n\n");
  /* Allocate memory for vectors */
  RHS=(double *) malloc(sizeof(double)*la);      /* Right-hand side vector */
  EX_SOL=(double *) malloc(sizeof(double)*la);   /* Analytical/exact solution */
  X=(double *) malloc(sizeof(double)*la);        /* Grid points */

  /* Initialize the problem: grid, RHS, and exact solution */
  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);                                /* Create uniform grid */
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);                     /* Set up RHS with BC */
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);  /* Compute exact solution */
  
  /* Write initial data to files for visualization */
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  /* Set up band storage parameters for tridiagonal matrix */
  kv=1;             /* Number of superdiagonals */
  ku=1;             /* Number of superdiagonals in original matrix */
  kl=1;             /* Number of subdiagonals */
  lab=kv+kl+ku+1;   /* Leading dimension of band storage */

  /* Allocate and initialize the coefficient matrix */
  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  printf("Solution with LAPACK\n");
  ipiv = (int *) calloc(la, sizeof(int));  /* Pivot indices for LU factorization */

  /* LU Factorization using LAPACK's general band factorization */
  /* DEBUG BEFORE FACTORIZATION */
  printf("DEBUG BEFORE: la=%d kl=%d ku=%d lab=%d NRHS=%d\n", la, kl, ku, lab, NRHS);
  for (int jj=0; jj < (la<6?la:6); ++jj){
    printf("DEBUG BEFORE AB col %d:", jj);
    for (int rr=0; rr < lab; ++rr) printf(" % .6e", AB[rr + jj*lab]);
    printf("\n");
  }

  clock_gettime(CLOCK_MONOTONIC, &t_start);

  if (IMPLEM == TRF) {
    dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }

  /* LU for tridiagonal matrix (can replace dgbtrf_) - custom implementation */
  if (IMPLEM == TRI) {
    dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  }


  /* DEBUG AFTER FACTORIZATION (ipiv and LU stored in AB) */
  printf("DEBUG AFTER: ipiv (first 10):");
  for (int kk=0; kk < (la<10?la:10); ++kk) printf(" %d", ipiv[kk]);
  printf("\n");
  for (int jj=0; jj < (la<6?la:6); ++jj){
    printf("DEBUG AFTER AB col %d:", jj);
    for (int rr=0; rr < lab; ++rr) printf(" % .6e", AB[rr + jj*lab]);
    printf("\n");
  }
  /* Back-substitution to solve the system after factorization */
  if (IMPLEM == TRI || IMPLEM == TRF){
    /* Solution (Triangular) - solve using the LU factors */
    if (info==0){
      dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
      if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
      printf("\n INFO = %d\n",info);
    }
  }

  /* Alternative: solve directly using dgbsv */
  if (IMPLEM == SV) {
    // TODO : use dgbsv
    printf("DEBUG: calling DGBSV with la=%d, kl=%d, ku=%d, lab=%d, NRHS=%d\n", la, kl, ku, lab, NRHS);
    dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    if (info != 0) {
      printf("\n INFO DGBSV (all-in-one solver) = %d\n", info);
    }
  }

  clock_gettime(CLOCK_MONOTONIC, &t_end);

  elapsed = (t_end.tv_sec - t_start.tv_sec)
        + 1e-9 * (t_end.tv_nsec - t_start.tv_nsec);

  printf("Elapsed time = %e s\n", elapsed);


  /* Write results to files */
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");  /* LU factors */
  write_xy(RHS, X, &la, "SOL.dat");  /* Solution at grid points (RHS now contains solution) */

  /* Relative forward error - compare numerical solution with exact solution */
  relres = relative_forward_error(RHS, EX_SOL, &la);
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  /* Free allocated memory */
  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");
}
