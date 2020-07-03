.SUFFIXES : .f .c .o

include ../makefile.in

INCLUDES = -I../INC -I../EXTERNAL/ITERSOL -I../EXTERNAL/DIRECTSOL
ifeq ($(INTEL), 1)
INCLUDES += -I${MKLROOT}/include
endif

# Object files
OBJS = cheblanNr.o lanvectors.o chebpol.o kpmdos.o lantrbnd.o landos.o lspolapprox.o pevsl.o \
       miscla.o parcsr.o parcsrmv.o parvec.o parvecs.o simpson.o spmat.o stats.o \
       spslicer.o utils.o vector.o \
       ../FORTRAN/pevsl_f90.o

LIB = ../libpevsl.a

# Rules
default: $(LIB)

%.f : %.o
	$(FORT) $(FLAGS) $(INCLUDES) -o $@ -c $<
%.o : %.c
	$(CC)   $(FLAGS) $(INCLUDES) -o $@ -c $<
$(LIB): $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f $(OBJS) $(LIB)
