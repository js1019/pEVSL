-----------------------------------------------------------------------
                  pEVSL (parallel EigenValue SLicing)
-----------------------------------------------------------------------
  Latest changes made on: 
-----------------------------------------------------------------------

main data structures for parallel computing

1. evsl_comm 
   This structure contains all process informations in the parallel environment
   and all the MPI communicators. There are 3 MPI_comm's that are 1) the global communicator
   (MPI_WORLD) 2), group communicator (for each slice), and 3) group lead communicator, i.e., 
   communicator for all the rank-0 processes in all the groups

2. evsl_parcsr 
   This structure is for parallel CSR matrices. Each process hold a block of rows, say Ai = A(r_i:r_{i+1}-1,:) 
   According to row/column partitioning, Ai is split into a local part (in the diagonal block) and an external part
   (off block diagonal), which are 2 CSR matrix, namely diag and offd.
   Currently, the only operation supported is parcsr_Matvec. A structure called commHandle is in evsl_parcsr, which takes care
   of all communications for matvec. Non-block communications are used in matvec to overlap computaions and communications.

3. evsl_parvec 
   This structure is for parallel vectors. 

FILES:
  * INC: [all header files]
    * pevsl.h               :   functino prototypes and constant definitions

    * pevsl_blaslapack.h    :   BLAS/LAPACK macros

    * pevsl_def.h           :   Macro definitions

    * pevsl_int.h           :   internal header with constant definitions and prototypes

    * pevsl_struct.h        :   pEVSL structs

  * SRC:	
    * cheblanNr.c	    :  Polynomial Filtered no-restart Lanczos

    * chebpol.c	    :  Computing and applying polynomial filters

    * kpmdos.c	    :  Compute DOS by KPM methods

    * landos.c	    :  Function to use Lanczos method for approximating DOS for the generalized eigenvalue problem

    * lantrbnd.c	    :  A more robust algorithm to give bounds of spectrum based on TR Lanczos

    * lspolapprox.c	:  Least squares polynomial approximation to a matrix function

    * miscla.c	    :  Miscellaneous la functions

    * parcsr.c	    :  Functinos related to parallel csr structure

    * parcsrmv.c	    :  Parallel csr matrix vector products

    * parvec.c	    :  Parallel vector related functions

    * parvecs.c	    :  Parallel multi-vector related functions

    * pevsl.c	        :  PEVSL interface functions

    * simpson.c	    :  Simpson integrater

    * spmat.c	        :  Space matrix operations

    * spslicer.c	    :  Spectrum slicing

    * stats.c	        :  Used to track various statistics (time taken by various operations)

    * utils.c	        :  Utility functions

    * vector.c	    :  Vector related functions 
- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

INSTALL:
     every effort is made so that the only thing to modify is the file
     makefile.in
          

      
      [0. need BLAS/LAPACK and MPI]
      1. modify makefile.in (may need to change MKLRoot, BLAS/LAPACK/MPI linking)
      2. make clean; make

      make in the directory SRC will create the library libpevsl.a

      In directory TESTS/ you will find a makefile to create 3 directorys
      * Gen_Lap              [Generalized Eigenvalue Laplacian]
      * MM                   [Matrix market form]
      * Lap                  [Laplacian]

Run:
```
      mpirun -np 10 ./test.ex -n 400 ./test.ex -nx 20 -ny 20 -nz 20 -nslices 5  -a 0.6 -b 1.2


      [5 slices, 2 processes for each slice]

= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 nx 20 ny 20 nz 20, nslices 5, a = 6.000000e-01 b = 1.200000e+00
= = = = = = = = = = = = = = = = SLICES  FOUND = = = = = = = = = = = = = = = = = = = =
 Slice 0: [6.000000000000000e-01, 7.431715857928962e-01]
 Slice 1: [7.431715857928962e-01, 8.719359679839922e-01]
 Slice 2: [8.719359679839922e-01, 9.904952476238122e-01]
 Slice 3: [9.904952476238122e-01, 1.098549274637319e+00]
 Slice 4: [1.098549274637319e+00, 1.200000000000000e+00]
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 Eigenvalue bounds: (0.000000e+00, 1.200000e+01)
 Eigenvalue count 126.44, estimated count for each slice 31
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
 Group  0 (NP=2): [6.00e-01, 7.43e-01], eval computed 18 (18), polydeg  61, its 210
 Group  1 (NP=2): [7.43e-01, 8.72e-01], eval computed 27 (27), polydeg  74, its 270
 Group  2 (NP=2): [8.72e-01, 9.90e-01], eval computed 27 (27), polydeg  86, its 270
 Group  3 (NP=2): [9.90e-01, 1.10e+00], eval computed 19 (19), polydeg  99, its 210
 Group  4 (NP=2): [1.10e+00, 1.20e+00], eval computed 27 (27), polydeg 110, its 270
= = = = = = = = = = = = = = = = = = = SLICE  0 = = = = = = = = = = = = = = = = = = = =
 Timing (sec):
   Create Comms  0.000628,   Create Parcsr 0.002328,   Create Slicer 0.000001
   Apply Slicer  0.245428,   FiltEig Solve 1.105189
 Memory (MB):
   Total 20.94,  Peak 18.60 
= = = = = = = = = = = = = = = = = = = SLICE  1 = = = = = = = = = = = = = = = = = = = =
 Timing (sec):
   Create Comms  0.000217,   Create Parcsr 0.002418,   Create Slicer 0.000001
   Apply Slicer  0.245796,   FiltEig Solve 1.689088
 Memory (MB):
   Total 21.52,  Peak 19.18 
= = = = = = = = = = = = = = = = = = = SLICE  2 = = = = = = = = = = = = = = = = = = = =
 Timing (sec):
   Create Comms  0.000205,   Create Parcsr 0.002302,   Create Slicer 0.000002
   Apply Slicer  0.245934,   FiltEig Solve 1.950552
 Memory (MB):
   Total 21.52,  Peak 19.18 
= = = = = = = = = = = = = = = = = = = SLICE  3 = = = = = = = = = = = = = = = = = = = =
 Timing (sec):
   Create Comms  0.000242,   Create Parcsr 0.002334,   Create Slicer 0.000002
   Apply Slicer  0.245882,   FiltEig Solve 1.717661
 Memory (MB):
   Total 21.01,  Peak 18.67 
= = = = = = = = = = = = = = = = = = = SLICE  4 = = = = = = = = = = = = = = = = = = = =
 Timing (sec):
   Create Comms  0.000216,   Create Parcsr 0.002417,   Create Slicer 0.000002
   Apply Slicer  0.245852,   FiltEig Solve 2.392627
 Memory (MB):
   Total 21.52,  Peak 19.18 
= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

      ```
