.SUFFIXES : .c .o

include ../../makefile.in

#FLAGS ?= -DUNIX -O3 -g -Wall 

INCLUDES = -I../../INC -I../COMMON -I../../EXTERNAL/DIRECTSOL

ifeq ($(INTEL), 1)
  OBJ_DIRSOL = ../../EXTERNAL/DIRECTSOL/pevsl_pardiso.o
else
  include ../../EXTERNAL/DIRECTSOL/makefile_mumps.in
  OBJ_DIRSOL = ../../EXTERNAL/DIRECTSOL/pevsl_mumps.o
  LIB_DIRSOL = $(LIB_MUMPS)
endif

# Object files
OBJS_CM = ../COMMON/comm.o ../COMMON/io.o ../COMMON/lapl.o
OBJS_PN  = LapPLanN.o $(OBJS_CM) $(OBJ_DIRSOL)

LIB = -L../../ -lpevsl 

ALLEXE = 
LIB_EXT = $(LIB_DIRSOL) -fopenmp $(LIBLAPACK) $(LIB0)

# intel Fortran lib
ifeq ($(INTEL), 1)
LIB_EXT += -lifcore
else
LIB_EXT += -lmpifort
endif

# Rules
default: LapPLanN.ex

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

# Exe
## polynomial
LapPLanN.ex: $(OBJS_PN) ../../libpevsl.a
	$(LINK) -o LapPLanN.ex $(OBJS_PN) $(LIB) $(LIB_EXT)

clean:
	rm -f *.o *.ex

cleanall:
	rm -f *.o *.ex OUT/*
