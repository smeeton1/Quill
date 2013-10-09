CC = gfortran
CFLAGS = 
LIB = -llapack 
# -lpacklib  -lpawlib -lpacklib -lmathlib -lgraflib -lgrafX11 \-lkernlib

OBJECTS = supperoperator.o          \
	  eigenSolve.o            \
	  phoenv.o               \
	  solver.o              \
	  until.o		\
	  special.o             \
	  utils.o               \
	  types.o		\
	  optimize.o		\
	  constants.o		\
	  amos.o		


%.o : src/%.f90
	$(CC) $(CFLAGS) -c $< $(LIB)

Photrun: $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LIB)


eigenSolve.o:

supperoperator.o:

amos.o:

types.o:

utils.o: types.o

constants.o: types.o

optimize.o: types.o utils.o

special.o: amos.o types.o utils.o constants.o optimize.o

solver.o: eigenSolve.o special.o

until.o: supperoperator.o

phoenv.o: until.o supperoperator.o eigenSolve.o solver.o 
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






clean:
	rm *o *mod Photrun
