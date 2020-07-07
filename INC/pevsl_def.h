#ifndef PEVSL_DEF_H
#define PEVSL_DEF_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

#define PEVSL_NOT_DEFINED -2
#define PEVSL_NOT_IMPLEMENT 11

#define PEVSL_MAX_LINE 2048

#define PEVSL_MAX(a, b) ((a) > (b) ? (a) : (b))
#define PEVSL_MIN(a, b) ((a) < (b) ? (a) : (b))

#define orthTol 1e-18

#define PI M_PI

#define PEVSL_ABORT(comm, errcode, msg) {\
    int pid; \
    MPI_Comm_rank(comm, &pid); \
    printf("PEVSL error (processor %d): %s \n", pid, msg);  \
    MPI_Abort(comm, errcode); \
}

#ifdef PEVSL_DEBUG
#define PEVSL_CHKERR(ierr) if (ierr) { PEVSL_ABORT(MPI_COMM_WORLD, 999, "Assertion failed") }
#else
#define PEVSL_CHKERR(ierr) 
#endif

/* memory management, alloc and free */

#define PEVSL_MALLOC(base, nmem, type) { \
  size_t nbytes = (nmem) * sizeof(type); \
  (base) = (type*) malloc(nbytes); \
  if ((base) == NULL) { \
    fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", nbytes); \
    fprintf(stdout, "Malloc at FILE %s, LINE %d, nmem %zu\n", __FILE__, __LINE__, (size_t) nmem); \
    MPI_Abort(MPI_COMM_WORLD, -1); \
  } \
}

#define PEVSL_CALLOC(base, nmem, type) { \
  size_t nbytes = (nmem) * sizeof(type); \
  (base) = (type*) calloc((nmem), sizeof(type)); \
  if ((base) == NULL) { \
    fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", nbytes); \
    fprintf(stdout, "Calloc at FILE %s, LINE %d, nmem %zu\n", __FILE__, __LINE__, (size_t) nmem); \
    MPI_Abort(MPI_COMM_WORLD, -1); \
  } \
}

#define PEVSL_REALLOC(base, nmem, type) {\
  size_t nbytes = (nmem) * sizeof(type); \
  (base) = (type*) realloc((base), nbytes); \
  if ((base) == NULL && nbytes > 0) { \
    fprintf(stdout, "EVSL Error: out of memory [%zu bytes asked]\n", nbytes); \
    fprintf(stdout, "Realloc at FILE %s, LINE %d, nmem %zu\n", __FILE__, __LINE__, (size_t) nmem); \
    MPI_Abort(MPI_COMM_WORLD, -1); \
  } \
}

#define PEVSL_FREE(base) {\
    if (base) { free(base); } \
    base = NULL; \
}

/*! These two MACROS are for running some code sequentially within an MPI_COMM
    PEVSL_SEQ_BEGIN(MPI_COMM);
    ...
    some code
    ...
    PEVSL_SEQ_END(MPI_COMM);
*/    
#define PEVSL_SEQ_BEGIN(MPI_COMM) { \
    int size, rank, __i; \
    MPI_Barrier(MPI_COMM); \
    MPI_Comm_size(MPI_COMM, &size); \
    MPI_Comm_rank(MPI_COMM, &rank); \
    /* sequential loop */ \
    for (__i=0; __i<size; __i++) { \
        if (__i == rank) { \

#define PEVSL_SEQ_END(MPI_COMM) \
        }\
        MPI_Barrier(MPI_COMM); \
    }\
}

/*!
  \def max(x,y)
  Computes the maximum of \a x and \a y.
*/
#define PEVSL_MAX(a, b) ((a) > (b) ? (a) : (b))

/*!
  \def min(x,y)
  Computes the minimum of \a x and \a y.
*/
#define PEVSL_MIN(a, b) ((a) < (b) ? (a) : (b))

/*! Fortran interface naming convention
 */
#define PEVSL_FORT(name) name ## _f90_

/*! max number of Gram–Schmidt process in orthogonalization
 */
#define NGS_MAX 2

#endif

