tool= mine

ifeq ($(tool), mine)
 CC = gfortran
 CFLAGS = 
 LIB = -llapack 
else
 CC = gfortran
 CFLAGS = 
 LIB =  -llapack -lf77blas -lcblas -latlas 
endif



# -lpacklib  -lpawlib -lpacklib -lmathlib -lgraflib -lgrafX11 \-lkernlib
 
include $(SLEPC_DIR)/conf/slepc_common 
PETSC_INCLUDE = -I$(PETSC_DIR)/include -I /usr/include/lam/
PETSC_ARCH_INCLUDE = -I$(PETSC_DIR)/$(PETSC_ARCH)/include



OBJECTS = supperoperator.o          \
	  eigenSolve.o            \
	  phoenv.o               \
	  solver.o              \
	  until.o		\
	  Densutil.o
# 	  special.o             \
# 	  utils.o               \
# 	  types.o		\
# 	  optimize.o		\
# 	  constants.o		\
# 	  amos.o		


%.o : src/%.F90
	$(CC) $(CFLAGS) -c $< $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)

Photrun: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)


eigenSolve.o:

supperoperator.o:

Densutil.o:

# amos.o:
# 
# types.o:
# 
# utils.o: types.o
# 
# constants.o: types.o
# 
# optimize.o: types.o utils.o
# 
# special.o: amos.o types.o utils.o constants.o optimize.o

until.o: supperoperator.o

solver.o: eigenSolve.o until.o

phoenv.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o
# all: Photrun
# 
# eigenSolve.o: eigenSolve.f90
# 	$(CC) $(CFLAGS) eigenSolve.f90 $(LIB)
# 
# supperoperator.o: supperoperator.f90
# 	$(CC) $(CFLAGS) supperoperator.f90
# 
# phoenv.o: phoenv.f90
# 	$(CC) $(CFLAGS) phoenv.f90 eigenSolve.o supperoperator.o  $(LIB)
# 
# Photrun: phoenv.o eigenSolve.o supperoperator.o
# 	$(CC) $(CFLAGS) phoenv.o eigenSolve.o supperoperator.o  $(LIB) -o Photrun






clean::
	rm *o *mod Photrun


.DEFAULT_GOAL := Photrun