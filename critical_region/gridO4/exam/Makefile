OBJ=FRG.f90 derivs.o dFdxO5.o chebft.o gapEq.o initial.o odeint.o phyphi.o rtbis.o spline.o splint.o  
F77=f95

exe: $(OBJ)
	$(F77) -o exe -O2 $(OBJ) $(LIBS)

derivs.o: derivs.f90
	$(F77)  -c -O2 derivs.f90

dFdxO5.o: dFdxO5.f90
	$(F77)  -c -O2 dFdxO5.f90

chebft.o: chebft.f90
	$(F77)  -c -O2 chebft.f90

gapEq.o: gapEq.f90
	$(F77)  -c -O2 gapEq.f90
	
initial.o: initial.f90
	$(F77)  -c -O2 initial.f90

odeint.o: odeint.f90
	$(F77)  -c -O2 odeint.f90

phyphi.o: phyphi.f90
	$(F77)  -c -O2 phyphi.f90
	
rtbis.o: rtbis.f90
	$(F77)  -c -O2 rtbis.f90

spline.o: spline.f90
	$(F77)  -c -O2 spline.f90

splint.o: splint.f90
	$(F77)  -c -O2 splint.f90

