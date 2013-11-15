tool= mine

ifeq ($(tool), mine)
 FLINKER = gfortran
 FLAGS = -w
 LIB = -llapack -I /usr/lib/openmpi/include/
else
 FLINKER = gfortran
 FLAGS = -w
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


All: Photrun 3edge

%.o : src/%.F90
	$(FLINKER) $(FLAGS) -c $< $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE)

Photrun: $(OBJECTS)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)
	
3edge: $(OBJECTS)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)


eigenSolve.o:

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

until.o: 

solver.o: eigenSolve.o until.o

supperoperator.o: until.o solver.o Densutil.o

phoenv.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o

3edge.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o


clean::
	rm *o *mod Photrun 3edge


.DEFAULT_GOAL := Photrun
