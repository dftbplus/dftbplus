FC=gfortran -O3 -fopenmp
OBJ =   functions.o magma_device_info.o
OUT=   libCountDevice.a
INCLUDES = -I/home/ccevallos/magma-2.3.0/include

$(OUT):  $(OBJ)
	ar rcs $(OUT) $(OBJ)	

%.o : %.f90  $(MODULES)
	$(FC) -c  $< 

%.o : %.cpp
	nvcc   $(INCLUDES) -c $@ $<

clean:
	rm -f *.mod *.o  *.a 
