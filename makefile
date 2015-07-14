MA?= mine
path?= obj/
BINDIR?= mod

ifeq ($(MA), mine)
 FLINKER = gfortran
 FLAGS = -w
 LIB =   -I /usr/lib/openmpi/include/ -I /usr/include -I /usr/local/lib/ -lfftw3 -llapack -lm
else
 FLINKER = ftn
 FLAGS = -w
 LIB =  -L$(FFTW_LIB) -lfftw3  
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
	  $(path)class.o          \
	  $(path)therm.o          \
	  $(path)phoenv.o               


OBJECTS3 = $(path)supperoperator.o          \
	   $(path)eigenSolve.o            \
	   $(path)solver.o              \
	   $(path)until.o		\
	   $(path)Densutil.o		\
	   $(path)func.o		\
	   $(path)class.o          \
	   $(path)directwalk.o          \
	   $(path)3edge.o		

OBJECTS4 = $(path)supperoperator.o          \
	   $(path)eigenSolve.o            \
	   $(path)solver.o              \
	   $(path)until.o		\
	   $(path)Densutil.o		\
	   $(path)func.o		\
	   $(path)class.o          \
	   $(path)directwalk.o          \
	   $(path)ranrun.o	
	   
OBJECTS5 = $(path)supperoperator.o          \
	   $(path)eigenSolve.o            \
	   $(path)solver.o              \
	   $(path)until.o		\
	   $(path)Densutil.o		\
	   $(path)func.o		\
	   $(path)class.o          \
	   $(path)directwalk.o          \
	   $(path)interact.o	

All: Photrun 3edge ranrun interact

$(path)%.o : src/%.F90
	$(FLINKER) $(FLAGS) -c $< $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) -o $@ -J$(BINDIR)
	
3edge: $(OBJECTS3)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS3) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)

ranrun: $(OBJECTS4)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS4) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)
	
Photrun: $(OBJECTS)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)
	
interact: $(OBJECTS5)
	$(FLINKER) $(FLAGS) -o $@ $(OBJECTS5) $(LIB) $(PETSC_INCLUDE) $(PETSC_ARCH_INCLUDE) $(SLEPC_INCLUDE) $(SLEPC_LIB)

eigenSolve.o:

Densutil.o:

until.o: 

func.o: eigenSolve.o until.o

solver.o: eigenSolve.o until.o

supperoperator.o: until.o solver.o Densutil.o func.o

class.o: until.o solver.o Densutil.o func.o eigenSolve.o

therm.o: until.o solver.o Densutil.o func.o supperoperator.o

Directwalk.o: until.o solver.o Densutil.o func.o supperoperator.o

phoenv.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o func.o class.o therm.o

3edge.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o  func.o class.o directwalk.o

ranrun.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o  func.o class.o directwalk.o

interact.o: until.o supperoperator.o eigenSolve.o solver.o Densutil.o  func.o class.o directwalk.o

clean::
	rm obj/*o mod/*.mod Photrun 3edge ranrun


.DEFAULT_GOAL := Photrun
