include ../makefile.in

ifeq ($(INTEL), 1)
  OBJ = DIRECTSOL/pevsl_pardiso.o ITERSOL/chebiter.o
  INCLUDES = -I../INC -I. -I$(MKLROOT)/include
else
  #include DIRECTSOL/makefile_mumps.in
  #OBJ = DIRECTSOL/pevsl_mumps.o ITERSOL/chebiter.o
  #INCLUDES = -I../INC -I. $(MUMPS_INC)
  #include DIRECTSOL/makefile_mumps.in
  OBJ = ITERSOL/chebiter.o
  INCLUDES = -I../INC -I. 
endif

%.o : %.c
	$(CC) $(FLAGS) $(INCLUDES) -o $@ -c $<

all : $(OBJ)

clean:
	rm -f DIRECTSOL/*.o ITERSOL/*.o

