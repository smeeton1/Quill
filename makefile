MA?= mine
path= obj/
BINDIR= mod

ifeq ($(MA), mine)
 FLINKER = gfortran
 FLAGS = -w
 LIB =   -I /usr/lib/openmpi/include/ -I /usr/include -I /usr/local/lib/ -lfftw3 -llapack -lm
else
 FLINKER = gfortran
 FLAGS = -w
 LIB =  -L/ivec/devel/intel/12.1/fftw/3.3.3/lib -lfftw3 -llapack -lf77blas -lcblas -latlas 
endif



# -lpacklib  -lpawlib -lpacklib -lmathlib -lgraflib -lgrafX11 \-lkernlib
 
include $(SLEPC_DIR)/conf/slepc_common 
PETSC_INCLUDE = -I$(PETSC_DIR)/include -I /usr/include/lam/
PETSC_ARCH_INCLUDE = -I$(PETSC_DIR)/$(PETSC_ARCH)/include



OBJECTS = $(path)eigenSolve.o            \
	  $(path)until.o		\
	  $(path)Densutil.o		\
	  $(path)solver.o              \
	  $(path)func.o			\
	  $(path)supperoperator.o          \
	  $(path)phoenv.o               
# 	  special.o             \
# 	  utils.o               \
# 	  types.o		\
# 	  optimize.o		\
# 	  constants.o		\
# 	  amos.o		

OBJECTS3 = $(path)supperoperator.o          \
	   $(path)eigenSolve.o            \
	   $(path)solver.o              \
	   $(path)until.o		\
	   $(path)Densutil.o		\
	   $(path)func.o		\
	   $(path)3edge.o		


All: Photrun 3edge

$(path)%.o : src/%.F90
	$(FLINKER) $(FLAGS) -c $< $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) -o $@ -J$(BINDIR)

Photrun: $(OBJECTS)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)
	
3edge: $(OBJECTS3)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS3) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)


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

func.o: eigenSolve.o until.o

solver.o: eigenSolve.o until.o

supperoperator.o: until.o solver.o Densutil.o func.o

phoenv.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o func.o

3edge.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o  func.o


clean::
	rm obj/*o mod/*.mod Photrun 3edge


.DEFAULT_GOAL := Photrun
