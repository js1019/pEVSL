#ifndef PEVSL_STRUCT_H
#define PEVSL_STRUCT_H

#include <mpi.h>
#include <stdlib.h>

/**
 * @brief sparse matrix format: the coordinate (COO) format, 0-based
 *
 * ir, jc, vv : triples for all nonzeros (of size nnz)
 */
typedef struct _pevsl_cooMat {
  int nrows, ncols, nnz, *ir, *jc;
  double *vv;
} pevsl_Coo;


/*! 
 * @brief sparse matrix format: the compressed sparse row (CSR) format, 0-based
 * 
 * 3-array variant: ia,ja,a, nnz == ia[nrows]
 */ 
typedef struct _pevsl_csrMat {
  int nrows, ncols, *ia, *ja;
  double  *a;
} pevsl_Csr;

/*!
 * @brief macro to obtain number of nonzeros of a CSR matrix 
 */
#define PEVSL_CSRNNZ(A) A->ia[A->nrows]

/*! 
 * @brief communication handle used in parallel matvec
 * 
 */ 
typedef struct _commHandle {
  int  num_proc_send_to;   /**< num of procs to send data to */
  int *proc_send_to;       /**< ranks of these procs */
  int *send_elmts_ids;     /**< local ids of all sending data */
  int *send_elmts_ptr;     /**< starting ptr of each piece of sending data, size of num_proc_send_to+1 */

  int  num_proc_recv_from; /**< num of procs to receive data from */
  int *proc_recv_from;     /**< ranks of these procs */
  int *recv_elmts_ptr;     /**< starting ptr of each piece of data reiceived [from each proc], size of num_proc_recv_from+1 */

  double *send_buf;
  double *recv_buf;
  
  MPI_Request *send_recv_requests; /**< communication request for MPI_Isend and Irecv */
  MPI_Status  *send_recv_status;   /**< status for MPI_Waitall */
} commHandle;


/*! 
 * @brief parallel CSR format
 * 
 */ 
typedef struct _pevsl_parcsr {
  /**< mpi communicator that this matrix resides on */
  MPI_Comm comm;

  /**< handle communications involved in parCSR matvec */
  commHandle *comm_handle;

  int nrow_global;
  int ncol_global;
  int nrow_local;
  int ncol_local;

  /* np: num of processors in comm, e.g., a group */
  /* array of length np+1, row_starts[i]: global index of the first row on proc i, row_starts[np] = nrow_global*/
  int *row_starts;
  /* array of length np+1, col_starts[i]: global index of the first column of the diag part on proc i, col_starts[np] = ncol_global */
  int *col_starts;

  /* my_first_row = row_starts[myrank],
   * my_first_col = col_starts[myrank]
   * my_row_range = [my_first_row : my_first_row + nrow_local-1]
   * my_col_range = [my_first_col : my_first_col + ncol_local-1]
   * 
   * For square matrix: row_range and col_range are often set to be the same
   * We reserve these for rectangular matrices (may not be useful pEVSL though)
   *
   * NOTE: if row_starts or col_starts == NULL, it means that the default 1D partitioning is used
   * Using routine pEVSL_Part1d can easily decide all the info regarding the partitioning,
   * see utils.c for details
   */

  int first_row;
  int first_col;

  /* diagonal    part of Ai    (local) : A(my_row_range, my_col_range) */
  pevsl_Csr *diag;

  /* off-diagonal part of Ai (external) : A(my_row_range, ~my_col_range) */
  pevsl_Csr *offd;

  /* global idx of the columns of offd, of size 'offd->ncols'
   * NOTE: columns are *local* indices. Use col_map_offd for mapping them back to global indices */
  int *col_map_offd;

} pevsl_Parcsr;


/*! 
 * @brief parallel vector
 */ 
typedef struct _pevsl_parvec {
  MPI_Comm comm; /**< MPI communicator */
  int n_global;  /**< number of elements in the global vector */
  int n_local;   /**< number of elements in the local portion */
  int n_first;   /**< the index of the first local element in the global vector
                      NOTE: this can be set as PEVSL_NOT_DEFINED */
  double *data;  /**< the data pointer */
  /* add for imagery parts JS 12/26/18 
  double *imag; */
} pevsl_Parvec;


/*! 
 * @brief parallel multiple vectors
 * global size: n_global x n_vecs
 *  local size: n_local  x n_vecs
 * in pEVSL, this is for storing Lanczos or Ritz vectors, so we often have n_global >> m
 * Thus, in the current code, we use only 1-D partitioning in rows for this data structure
 */
typedef struct _pevsl_parvecs {
  MPI_Comm comm; /**< MPI communicator */
  int n_global;  /**< number of elements in the global vector */
  int n_local;   /**< number of elements in the local portion */
  int n_first;   /**< the index of the first local element in the global vector
                      NOTE: this can be set as PEVSL_NOT_DEFINED */
  int n_vecs;    /**< the number of vectors */
  int ld;        /**< leading dimension (must be >= n_local) */
  double *data;  /**< the data pointer */
} pevsl_Parvecs;


/*!
 * @brief  parameters for polynomial filter 
 *
 * default values are set by set_pol_def
 */
typedef struct _pevsl_polparams {
  /** @name input to find_pol */
  /**@{*/
  int max_deg;        /**< max allowed degree */
  int min_deg ;       /**< min allowed degree */
  int damping;        /**< 0 = no damping, 1 = Jackson, 2 = Lanczos */
  double thresh_ext;  /**< threshold for accepting polynom. for end intervals */
  double thresh_int;  /**< threshold for interior intervals */
  double intvtol;     /**< cut-off point of middle interval */
  /**@}*/

  /** @name output from find_pol */
  /**@{*/
  int    type;       /**< type of the filter: 
                          0: middle interval, 1: left interval, 2: right interval */
  double *mu;         /**< coefficients. allocation done by set_pol */
  double cc;          /**< center of interval - used by chebAv */
  double dd;          /**< half-width of interval - used by chebAv */
  double gam;         /**< center of delta function used */
  double bar;         /**< p(theta)>=bar indicates a wanted Ritz value */
  /**@}*/

  /** @name both input to and output from find_pol */
  /**@{*/
  int deg ;           /**< if deg == 0 before calling find_deg then
                        the polynomial degree is  computed
                        internally. Otherwise it is of degree deg.
                        [and  thresh_ext and thresh_int are not used]
                        default value=0, set by call to set_pol_def */
  /**@}*/
} pevsl_Polparams;


/**
 * @brief linear solver function prototype: [complex version]
 * which is used for solving system with A-SIGMA B 
 * n  is the size  of the system,  br, bz are  the right-hand
 * side (real and  imaginary parts of complex vector),  xr, xz will
 * be the  solution (complex vector),  and "data" contains  all the
 * data  needed  by  the  solver. 
 */
/*
typedef void (*SVFuncC)(pevsl_Parvec *br, pevsl_Parvec *bz, 
                        pevsl_Parvec *xr, pevsl_Parvec *xz, 
                        void *data);
*/
typedef void (*SVFuncC)(double *br, double *bz, 
                        double *xr, double *xz, 
                        void *data);

/** 
 * @brief function prototype for applying the solve B x = b 
 * [the most general form]
 */
/* 
typedef void (*SVFunc)(pevsl_Parvec *b, pevsl_Parvec *x, void *data);
*/
typedef void (*SVFunc)(double *b, double *x, void *data);

typedef void (*ZSVFunc)(double *br, double *bi, double *xr, double *xi, void *data);

/**
 * @brief matvec function prototype 
 */
/*
typedef void (*MVFunc)(pevsl_Parvec *x, pevsl_Parvec *y, void *data);
*/
typedef void (*MVFunc)(double *x, double *y, void *data);

typedef void (*ZMVFunc)(double *xr, double *xi, double *yr, double *yi, void *data);

/*!
 * @brief user-provided Mat-Vec function and data for y = A * x or y = B * x
 *
 */
typedef struct _pevsl_Matvec {
  MVFunc func;         /**< function pointer */
  void *data;          /**< data */
} pevsl_Matvec;

/*!
 * @brief user-provided function and data for solving B x = b
 *
 */
typedef struct _pevsl_Bsol {
  SVFunc func;       /**< function pointer */
  void *data;        /**< data */
} pevsl_Bsol;

/*!
 * @brief user-provided function for solving L^{T} x = b
 *
 */
typedef struct _pevsl_LtSol {
  SVFunc func;       /**< function pointer */
  void *data;        /**< data */
} pevsl_Ltsol;

/*! JS 01/02/2019 add function for complex solve */
/*!
 * @brief user-provided Mat-Vec function and data for y = A * x or y = B * x
 */
typedef struct _pevsl_ZMatvec {
  ZMVFunc func;         /**< function pointer */
  void *data;          /**< data */
} pevsl_ZMatvec;
/*!
 * @brief user-provided function and data for solving complex B x = b
 */
typedef struct _pevsl_ZBsol {
  ZSVFunc func;       /**< function pointer */
  void *data;        /**< data */
} pevsl_ZBsol;



/*!
 * @brief timing and memory statistics of pEVSL
 *
 */
typedef struct _pevsl_Stat {
  /* timing [level-1 funcs] */
  double t_setBsv;
  double t_setASigBsv;
  double t_eigbounds;
  double t_dos;
  double t_iter;
  /* timing [level-2 funcs] */
  double t_mvA;
  double t_mvB;
  double t_svB;
  double t_svLT;
  double t_svASigB;
  double t_reorth;
  double t_eig;
  double t_blas;
  double t_ritz;
  double t_polAv;
  double t_ratAv;
  double t_sth;
  size_t n_mvA;
  size_t n_mvB;
  size_t n_svB;
  size_t n_svLT;
  size_t n_svASigB;
  size_t n_polAv;
  size_t n_ratAv;
  /* memory */
  size_t alloced;
  size_t alloced_total;
  size_t alloced_max;
} pevsl_Stat;

/*!
 * @brief wrapper of all global variables in pEVSL
 *
 */
typedef struct _pevsl_Data {
  /** We keep some useful information in this struct
   *  N is the global size of the problem, i.e, size of A (or B)
   *  n is the local size, number of rows owned by this MPI rank.
   *  If the matrix is assumed to be partitioned into blocks of consecutive rows,
   *  nfirst is the first row owned.
   *  NOTE that currently Parcsr matrix is always assumed to have such partitioning. 
   *  More general partitionings are NOT yet supported
   *  BUT, we leave the option to have such non-consecutive row partitioning. 
   *  The users can implement their own Amv, Bmv, Bsol routines with arbitrary partitionings,
   *  and only set N and n but leave nfirst as PEVSL_NOT_DEFINED
   */
  MPI_Comm      comm;       /**< MPI communicator where this instance of pEVSL lives;
                                 all parallel matrices and vectors should the SAME MPI_Comm */
  int           N;          /**< global size of matrix A and B */
  int           n;          /**< local size of matrix A and B */
  int           nfirst;     /**< the first local row and column */
  int           ifGenEv;    /**< if it is a generalized eigenvalue problem */
  pevsl_Matvec *Amv;        /**< external matvec routine and the associated data for A */
  pevsl_Matvec *Bmv;        /**< external matvec routine and the associated data for B */
  pevsl_Bsol   *Bsol;       /**< external function and data for B solve */
  pevsl_Ltsol  *LTsol;      /**< external function and data for LT solve */
  pevsl_Stat   *stats;      /**< timing and memory statistics of pEVSL */


  /* JS 01/02/19 add additional functions for complex systems*/
  pevsl_ZMatvec *ZAmv;  
  pevsl_ZMatvec *ZBmv; 
  pevsl_ZBsol   *ZBsol; 


  int            nev_computed;    /**< Used in Fortran interface:
                                       hold the points of last computed results */
  double        *eval_computed;
  pevsl_Parvecs *evec_computed;
  /* add for imagery parts JS 12/26/18 */ 
  pevsl_Parvecs *evec_imag_computed;
 
  double        sigma_mult; /** multiplier for sigma in LanDOS*/
} pevsl_Data;

#endif
