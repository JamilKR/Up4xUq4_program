FC = gfortran
FFLAGS = -Wall -O3

test_src = MOD_matfun.f90 MOD_Up4.f90 test.f90
test_obj = $(addsuffix .o, $(basename ${test_src}))

%.o: %.f90 $(test_src)
	$(FC)  $(FFLAGS) -c -o $@ $<

test: $(test_obj)
	$(FC) $(FFLAGS) -o $@".out" $^

clean:
	@rm *.o *.mod *.out
