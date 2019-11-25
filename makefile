
FC=mpif90
FFLAGS=-I$(STELLOPT_PATH)/LIBSTELL/Release/
LDFLAGS=$(STELLOPT_PATH)/LIBSTELL/Release/libstell.a

OBJ=spline_test.o

%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)

xspline_test: $(OBJ)
	$(FC) -o $@ $^ $(LDFLAGS) $(LIBS)

run: xspline_test
	@mpirun -np 1 xspline_test >& log_Fspl.txt
