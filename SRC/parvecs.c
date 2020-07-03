#include "pevsl_int.h"
/**
 * @file parvecs.c
 * @brief Parallel multi-vector related functions
 */

/*!
 * @brief Create a parallel multi-vector struct without allocating memory for data
 * @param[in] nglobal Number of global elements
 * @param[in] nvecs number of vecs
 * @param[in] ld lading dimension
 * @param[in] nlocal Number local elements
 * @param[in] nfirst Index of first local element
 * @param[in] comm Communicator
 * @param[out] x Vectors
 * @param[in] data Data
 */
void pEVSL_ParvecsCreateShell(int nglobal, int nvecs, int ld, int nlocal, int nfirst,
                              MPI_Comm comm, pevsl_Parvecs *x, double *data) {
  x->comm = comm;
  x->n_global = nglobal;
  x->n_local = nlocal;
  x->n_first = nfirst < 0 ? PEVSL_NOT_DEFINED : nfirst;
  x->n_vecs = nvecs;
  x->ld = ld;
  PEVSL_CHKERR(ld < nlocal);
  x->data = data;
}

/*!
 * @brief Create a parallel multi-vector struct
 * @param[in] nglobal Number of global elements
 * @param[in] nvecs Number of vectors
 * @param[in] ld leading dimension
 * @param[in] nlocal Number local elements
 * @param[in] nfirst Index of first local element
 * @param[in] comm Communicator
 * @param[out] x Vectors
 */
void pEVSL_ParvecsCreate(int nglobal, int nvecs, int ld, int nlocal, int nfirst,
                         MPI_Comm comm, pevsl_Parvecs *x) {

  pEVSL_ParvecsCreateShell(nglobal, nvecs, ld, nlocal, nfirst, comm, x, NULL);

  /* to prevent integer overflow, cast it to size_t before multiplication */
  size_t alloc = ((size_t) ld) * ((size_t) nvecs);
  PEVSL_MALLOC(x->data, alloc, double);
}

/*!
 * @brief Creates multi-vectors of the same type as an existing vector
 * @param[in] x input vector
 * @param[in] nvecs Number of vectors
 * @param[in] ld leading dimension
 * @param[in] y output vectors
 */
void pEVSL_ParvecsDuplParvec(pevsl_Parvec *x, int nvecs, int ld, pevsl_Parvecs *y) {

  pEVSL_ParvecsCreate(x->n_global, nvecs, ld, x->n_local, x->n_first, x->comm, y);
}

/*!
 * @brief Destroy a Parvecs struct
 */
void pEVSL_ParvecsFree(pevsl_Parvecs *x) {

  PEVSL_FREE(x->data);
}

/*!
 * @brief enables access vector i of a Parvecs as a Parvec
 * @param[in] X input vectors
 * @param[in] i Requested element
 * @param[out] x The ith vector of X
 */
void pEVSL_ParvecsGetParvecShell(pevsl_Parvecs *X, int i, pevsl_Parvec *x) {

  PEVSL_CHKERR(i >= X->n_vecs);
  double *data = i < 0 ? NULL : X->data + ((size_t) X->ld) * ((size_t) i);
  pEVSL_ParvecCreateShell(X->n_global, X->n_local, X->n_first, X->comm, x, data);
}

/*!
 * @brief Perform GEMV with a Parvecs and an double array to get a Parvec
 * y = alp * A(:,0:nvecs-1) * x + bet * y
 * @param[in] alp alp
 * @param[in] A input vectors
 * @param[in] nvecs Number of vectors
 * @param[in] x input vector (double array)
 * @param[in] bet bet
 * @param[out] y Output vector
 */
void pEVSL_ParvecsGemv(double alp, pevsl_Parvecs *A, int nvecs, double *x,
                       double bet, pevsl_Parvec *y) {
  PEVSL_CHKERR(A->n_global != y->n_global);
  PEVSL_CHKERR(A->n_local != y->n_local);
  PEVSL_CHKERR(A->n_first != y->n_first);
  char cN = 'N';
  int one = 1;
  DGEMV(&cN, &A->n_local, &nvecs, &alp, A->data, &A->ld, x, &one, &bet,
        y->data, &one);
}


/*!
 * @brief Perform GEMV with the transpose of a Parvecs and a Parvec
 * y = alp * A(:,0:nvecs-1)^T * x + bet * y
 * y is the same array stored on every proc
 * w is the workspace of size 2*nvecs
 * @param[in] alp alp
 * @param[in] A A
 * @param[in] nvecs Number of vectors
 * @param[in] x input vector
 * @param[in] bet bet
 * @param[in, out] y
 *
 * @param[in, out] w work vector
 */
void pEVSL_ParvecsGemtvWithWspace(double alp, pevsl_Parvecs *A, int nvecs, pevsl_Parvec *x,
                                  double bet, double *y, double *w) {
  PEVSL_CHKERR(A->n_global != x->n_global);
  PEVSL_CHKERR(A->n_local != x->n_local);
  PEVSL_CHKERR(A->n_first != x->n_first);
  char cT = 'T';
  int one = 1;
  double done = 1.0, dzero = 0.0, *w1, *w2;

  w1 = w;
  w2 = w + nvecs;

  DGEMV(&cT, &A->n_local, &nvecs, &done, A->data, &A->ld, x->data, &one, &dzero,
        w1, &one);

  MPI_Allreduce(w1, w2, nvecs, MPI_DOUBLE, MPI_SUM, A->comm);

  DSCAL(&nvecs, &bet, y, &one);
  DAXPY(&nvecs, &alp, w2, &one, y, &one);
}


/*!
 * @brief Perform GEMV with the transpose of a Parvecs and a Parvec
 * y = alp * A(:,0:nvecs-1)^T * x + bet * y
 * y is the same array stored on every proc
 * w is the workspace of size 2*nvecs
 * @param[in] alp alp
 * @param[in] A A
 * @param[in] nvecs Number of vectors
 * @param[in] x input vector
 * @param[in] bet bet
 * @param[in, out] y
 */
void pEVSL_ParvecsGemtv(double alp, pevsl_Parvecs *A, int nvecs, pevsl_Parvec *x,
                        double bet, double *y) {
  PEVSL_CHKERR(A->n_global != x->n_global);
  PEVSL_CHKERR(A->n_local != x->n_local);
  PEVSL_CHKERR(A->n_first != x->n_first);
  char cT = 'T';
  int one = 1;
  double done = 1.0, dzero = 0.0, *w1, *w2;

  PEVSL_MALLOC(w1, 2*nvecs, double);
  w2 = w1 + nvecs;

  DGEMV(&cT, &A->n_local, &nvecs, &done, A->data, &A->ld, x->data, &one, &dzero,
        w1, &one);

  MPI_Allreduce(w1, w2, nvecs, MPI_DOUBLE, MPI_SUM, A->comm);

  DSCAL(&nvecs, &bet, y, &one);
  DAXPY(&nvecs, &alp, w2, &one, y, &one);

  PEVSL_FREE(w1);
}

