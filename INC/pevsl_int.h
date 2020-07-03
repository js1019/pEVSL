/*
  This file is the internal header file of pEVSL that contains data structures, 
  constant definitions, and internal function prototypes 
*/

#ifndef PEVSL_INT_H
#define PEVSL_INT_H

#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <float.h>

#include "pevsl_struct.h"
#include "pevsl_def.h"
#include "pevsl_blaslapack.h"
#include "pevsl.h"

/* chebpol.c */
int pEVSL_ChebAv(pevsl_Data *pevsl, pevsl_Polparams *pol, pevsl_Parvec *v, pevsl_Parvec *y, pevsl_Parvec *w);
// JS 020619 
int pEVSL_ZChebAv(pevsl_Data *pevsl, pevsl_Polparams *pol, 
                  pevsl_Parvec *vr, pevsl_Parvec *vi, 
                  pevsl_Parvec *yr, pevsl_Parvec *yi, 
                  pevsl_Parvec *wr, pevsl_Parvec *wi);

int dampcf(int m, int damping, double *jac);
int chebxPltd(int m, double *mu, int npts, double *xi, double *yi);

/* miscla.c */
int SymmTridEig(pevsl_Data *pevsl, double *eigVal, double *eigVec, int n, const double *diag, const double *sdiag);

int SymmTridEigS(pevsl_Data *pevsl, double *eigVal, double *eigVec, int n, double vl, double vu, int *nevO, const double *diag, const double *sdiag);

void SymEigenSolver(pevsl_Data *pevsl, int n, double *A, int lda, double *Q, int ldq, double *lam);

void CGS_DGKS(pevsl_Data *pevsl, int k, int i_max, pevsl_Parvecs *Q, pevsl_Parvec *v, double *nrmv, double *w);

void CGS_DGKS2(pevsl_Data *pevsl, int k, int i_max, pevsl_Parvecs *V, pevsl_Parvecs *Z, pevsl_Parvec *v, double *w);

/* add JS */
void CGS_ZDGKS(pevsl_Data *pevsl, int k, int i_max, 
              pevsl_Parvecs *Qr, pevsl_Parvecs *Qi, 
              pevsl_Parvec  *vr, pevsl_Parvec  *vi,
              double *nrmv, double *wr, double *wi);

void CGS_ZDGKS2(pevsl_Data *pevsl, int k, int i_max, 
         pevsl_Parvecs *Vr, pevsl_Parvecs *Vi, pevsl_Parvecs *Zr, pevsl_Parvecs *Zi, 
         pevsl_Parvec  *vr, pevsl_Parvec  *vi, double *wr, double *wi);

/* parcsrmv.c */
void pEVSL_ParcsrMatvec0(double *x, double *y, void *data);

/* simpson.c */
void simpson(double* xi, double* yi, int npts);

/* utils.c */
int pEVSL_BinarySearchInterval(int *x, int n, int key);

/*------------------- inline functions */
/** 
 * @brief fprintf, only for rank == 0
 * 
 * */
static inline void pEVSL_fprintf0(int rank, FILE *fp, const char *format, ...) {

  if (rank) {
    return;
  }
  va_list args;
  va_start(args, format);
  vfprintf(fp, format, args);
  va_end(args);
}

/** 
 * @brief Perform matrix-vector product y = A * x
 * 
 * */
static inline void pEVSL_MatvecA(pevsl_Data   *pevsl_data, 
                                 pevsl_Parvec *x, 
                                 pevsl_Parvec *y) {

  PEVSL_CHKERR(!pevsl_data->Amv);
     
  PEVSL_CHKERR(pevsl_data->N != x->n_global);
  PEVSL_CHKERR(pevsl_data->n != x->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != x->n_first);
  PEVSL_CHKERR(pevsl_data->N != y->n_global);
  PEVSL_CHKERR(pevsl_data->n != y->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != y->n_first);

  double tms = pEVSL_Wtime();

  pevsl_data->Amv->func(x->data, y->data, pevsl_data->Amv->data);
  
  double tme = pEVSL_Wtime();
  pevsl_data->stats->t_mvA += tme - tms;
  pevsl_data->stats->n_mvA ++;
}


/** 
 * @brief Perform matrix-vector product y = B * x
 * 
 * */
static inline void pEVSL_MatvecB(pevsl_Data   *pevsl_data,
                                 pevsl_Parvec *x, 
                                 pevsl_Parvec *y) {

  PEVSL_CHKERR(!pevsl_data->Bmv);
  
  PEVSL_CHKERR(pevsl_data->N != x->n_global);
  PEVSL_CHKERR(pevsl_data->n != x->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != x->n_first);
  PEVSL_CHKERR(pevsl_data->N != y->n_global);
  PEVSL_CHKERR(pevsl_data->n != y->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != y->n_first);

  double tms = pEVSL_Wtime();
  
  pevsl_data->Bmv->func(x->data, y->data, pevsl_data->Bmv->data);
  
  double tme = pEVSL_Wtime();
  pevsl_data->stats->t_mvB += tme - tms;
  pevsl_data->stats->n_mvB ++;
}

/**
* @brief y = B \ x
* This is the solve function for the matrix B in pevsl_Data
*/
static inline void pEVSL_SolveB(pevsl_Data   *pevsl_data,
                                pevsl_Parvec *x, 
                                pevsl_Parvec *y) {

  PEVSL_CHKERR(!pevsl_data->Bsol);
 
  PEVSL_CHKERR(pevsl_data->N != x->n_global);
  PEVSL_CHKERR(pevsl_data->n != x->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != x->n_first);
  PEVSL_CHKERR(pevsl_data->N != y->n_global);
  PEVSL_CHKERR(pevsl_data->n != y->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != y->n_first);

  double tms = pEVSL_Wtime();
  
  pevsl_data->Bsol->func(x->data, y->data, pevsl_data->Bsol->data);
  
  double tme = pEVSL_Wtime();
  pevsl_data->stats->t_svB += tme - tms;
  pevsl_data->stats->n_svB ++;
}

/* JS introduce complex version */
/** 
 * @brief Perform complex matrix-vector product y = A * x
 * 
 * */
static inline void pEVSL_ZMatvecA(pevsl_Data   *pevsl_data, 
                 pevsl_Parvec *xr, pevsl_Parvec *xi, 
                 pevsl_Parvec *yr, pevsl_Parvec *yi) {

  PEVSL_CHKERR(!pevsl_data->ZAmv);
     
  PEVSL_CHKERR(pevsl_data->N != xr->n_global);
  PEVSL_CHKERR(pevsl_data->n != xr->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != xr->n_first);
  PEVSL_CHKERR(pevsl_data->N != yr->n_global);
  PEVSL_CHKERR(pevsl_data->n != yr->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != yr->n_first);

  double tms = pEVSL_Wtime();

  pevsl_data->ZAmv->func(xr->data, xi->data,
                         yr->data, yi->data, pevsl_data->ZAmv->data);
  
  double tme = pEVSL_Wtime();
  pevsl_data->stats->t_mvA += tme - tms;
  pevsl_data->stats->n_mvA ++;
}

/** 
 * @brief Perform complex matrix-vector product y = B * x
 * 
 * */
static inline void pEVSL_ZMatvecB(pevsl_Data   *pevsl_data,
                       pevsl_Parvec *xr, pevsl_Parvec *xi, 
                       pevsl_Parvec *yr, pevsl_Parvec *yi) {

  PEVSL_CHKERR(!pevsl_data->ZBmv);
  
  PEVSL_CHKERR(pevsl_data->N != xr->n_global);
  PEVSL_CHKERR(pevsl_data->n != xr->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != xr->n_first);
  PEVSL_CHKERR(pevsl_data->N != yr->n_global);
  PEVSL_CHKERR(pevsl_data->n != yr->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != yr->n_first);

  double tms = pEVSL_Wtime();
  
  pevsl_data->ZBmv->func(xr->data, xi->data,
                        yr->data, yi->data,pevsl_data->ZBmv->data);
  
  double tme = pEVSL_Wtime();
  pevsl_data->stats->t_mvB += tme - tms;
  pevsl_data->stats->n_mvB ++;
}

/**
* @brief complex y = B \ x
* This is the solve function for the matrix B in pevsl_Data
*/
static inline void pEVSL_ZSolveB(pevsl_Data   *pevsl_data,
                       pevsl_Parvec *xr, pevsl_Parvec *xi, 
                       pevsl_Parvec *yr, pevsl_Parvec *yi) {

  PEVSL_CHKERR(!pevsl_data->ZBsol);
 
  PEVSL_CHKERR(pevsl_data->N != xr->n_global);
  PEVSL_CHKERR(pevsl_data->n != xr->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != xr->n_first);
  PEVSL_CHKERR(pevsl_data->N != yi->n_global);
  PEVSL_CHKERR(pevsl_data->n != yi->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != yi->n_first);

  double tms = pEVSL_Wtime();
  
  pevsl_data->ZBsol->func(xr->data, xi->data, yr->data, yi->data, pevsl_data->ZBsol->data);
  
  double tme = pEVSL_Wtime();
  pevsl_data->stats->t_svB += tme - tms;
  pevsl_data->stats->n_svB ++;
}



/**
* @brief y = LT \ x or y = SQRT(B) \ x
* This is the solve function for the matrix B in pevsl_Data
*/
static inline void pEVSL_SolveLT(pevsl_Data   *pevsl_data,
                                 pevsl_Parvec *x, 
                                 pevsl_Parvec *y) {

  PEVSL_CHKERR(!pevsl_data->LTsol);
 
  PEVSL_CHKERR(pevsl_data->N != x->n_global);
  PEVSL_CHKERR(pevsl_data->n != x->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != x->n_first);
  PEVSL_CHKERR(pevsl_data->N != y->n_global);
  PEVSL_CHKERR(pevsl_data->n != y->n_local);
  PEVSL_CHKERR(pevsl_data->nfirst != y->n_first);

  double tms = pEVSL_Wtime();
  
  pevsl_data->LTsol->func(x->data, y->data, pevsl_data->LTsol->data);
  
  double tme = pEVSL_Wtime();
  pevsl_data->stats->t_svLT += tme - tms;
  pevsl_data->stats->n_svLT ++;
}

/**
 * @brief check if an interval is valid 
 **/
static inline int check_intv(double *intv, FILE *fstats) {
  /* intv[4]: ( intv[0], intv[1] ) is the inteval of interest
   *          ( intv[2], intv[3] ) is the spectrum bounds
   * return   0: ok
   *        < 0: interval is invalid
   */
  double a=intv[0], b=intv[1], lmin=intv[2], lmax=intv[3];
  if (a >= b) {
    fprintf(fstats, " error: invalid interval (%e, %e)\n", a, b);
    return -1;
  }
  
  if (a >= lmax || b <= lmin) {
    fprintf(fstats, " error: interval (%e, %e) is outside (%e %e) \n", 
            a, b, lmin, lmax);
    return -2;
  } 
  
  return 0;
}

#endif

