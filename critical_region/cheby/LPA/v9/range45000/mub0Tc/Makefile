OBJ=FRG.f90 derivs.o chebft.o fun.o nfdx.o initial.o odeint.o phyphi.o threshold_fun.o 
F77=gfortran

exe: $(OBJ)
	$(F77) -o exe -O2 $(OBJ) $(LIBS)

derivs.o: derivs.f90
	$(F77)  -c -O2 derivs.f90

chebft.o: chebft.f90
	$(F77)  -c -O2 chebft.f90

fun.o: fun.f90
	$(F77)  -c -O2 fun.f90

nfdx.o: nfdx.f90
	$(F77)  -c -O2 nfdx.f90

initial.o: initial.f90
	$(F77)  -c -O2 initial.f90

odeint.o: odeint.f90
	$(F77)  -c -O2 odeint.f90

phyphi.o: phyphi.f90
	$(F77)  -c -O2 phyphi.f90

threshold_fun.o: threshold_fun.f90
	$(F77)  -c -O2 threshold_fun.f90
