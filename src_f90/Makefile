################################################################################
DISLIN=${HOME}/dislin
f90comp=gfortran -Wall -O2
lib_flag = -L${DISLIN} -ldislin_d
# lib_flag = ${DISLIN}/libdislin_d.a
inc_flag = -I${DISLIN}/gf/real64
################################################################################
main.out : main.f90 mydata.o ode.o intersection.o
	$(f90comp) -o main.out main.f90 mydata.o ode.o intersection.o $(lib_flag)
mydata.o : mydata.f90
	$(f90comp) -c mydata.f90
# output.o : output.f90 mydata.f90
# 	$(f90comp) -c output.f90 mydata.f90
# interpolation.o : interpolation.f90 mydata.f90
# 	$(f90comp) -c interpolation.f90 mydata.f90
ode.o : ode.f90 mydata.f90
	$(f90comp) -c ode.f90 mydata.f90
# root_find.o : root_find.f90 mydata.f90 ode.f90
# 	$(f90comp) -c root_find.f90 mydata.f90 ode.f90
intersection.o : intersection.f90 mydata.f90
	$(f90comp) $(inc_flag) -c intersection.f90 mydata.f90
################################################################################
clean:
	rm *.mod *.o main.out
################################################################################
run:
	LD_LIBRARY_PATH=${DISLIN}:${LD_LIBRARY_PATH} && ./main.out
