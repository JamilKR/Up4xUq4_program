FC = gfortran
FFLAGS = -g -Wall -O3 -fcheck=all

MODLIB = -I/usr/local/lib/lapack95/lapack95_modules
OPTLIB = -lblas -llapack /usr/local/lib/lapack95/lapack95.a

test_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 test.f90
test_obj = $(addsuffix .o, $(basename ${test_src}))

test_2_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 test_2.f90
test_2_obj = $(addsuffix .o, $(basename ${test_2_src}))

%.o: %.f90 $(test_src) $(test_2_src)
	$(FC) $(FFLAGS) $(MODLIB) -c -o $@ $<

test: $(test_obj)
	$(FC) $(FFLAGS) $(MODLIB) -o $@".out" $^ $(OPTLIB)

test_2: $(test_2_obj)
	$(FC) $(FFLAGS) $(MODLIB) -o $@".out" $^ $(OPTLIB)

clean:
	@rm *.o *.mod *.out
