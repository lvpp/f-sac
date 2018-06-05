f-sac: Activity.o CompsData.o Allocation.o RaoultPressure.o main.o
	gfortran -o f-sac Activity.o CompsData.o Allocation.o RaoultPressure.o main.o
variables.mod: Allocation.o Allocation.f90
	gfortran -c Allocation.f90
deltaw_hb_data.mod: Allocation.o Allocation.f90
	gfortran -c Allocation.f90
Allocation.o: Allocation.f90
	gfortran -c Allocation.f90
Activity.o: variables.mod Activity.f90
	gfortran -c Activity.f90
CompsData.o: variables.mod CompsData.f90
	gfortran -c CompsData.f90
RaoultPressure.o: variables.mod RaoultPressure.f90
	gfortran -c RaoultPressure.f90
main.o: variables.mod main.f90
	gfortran -c main.f90
clean:
	rm -f f-sac variables.mod deltaw_hb_data.mod Activity.o CompsData.o Allocation.o RaoultPressure.o main.o
