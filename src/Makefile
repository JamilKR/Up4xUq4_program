FC = gfortran
FFLAGS = -g -Wall -O3 -fcheck=all

MODLIB = -I/usr/local/lib/lapack95/lapack95_modules
OPTLIB = -lblas -llapack /usr/local/lib/lapack95/lapack95.a

test_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 test.f90
test_obj = $(addsuffix .o, $(basename ${test_src}))

test_2_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 test_2.f90
test_2_obj = $(addsuffix .o, $(basename ${test_2_src}))

test_3_src = MOD_matfun.f90  MOD_Up_x_Uq.f90 test_3.f90
test_3_obj = $(addsuffix .o, $(basename ${test_3_src}))

%.o: %.f90 $(test_src) $(test_2_src) $(test_3_src)
	$(FC) $(FFLAGS) $(MODLIB) -c -o $@ $<

all: help test test_2 test_3

help:
	@echo
	@echo   "\e[42mDemasiadas expectativas...!\e[0m"
	@echo
	@echo   "\e[44mToo many expectations...!\e[0m"
	@echo

test: $(test_obj)
	$(FC) $(FFLAGS) $(MODLIB) -o $@".out" $^ $(OPTLIB)

test_2: $(test_2_obj)
	$(FC) $(FFLAGS) $(MODLIB) -o $@".out" $^ $(OPTLIB)

test_3: $(test_3_obj)
	$(FC) $(FFLAGS) $(MODLIB) -o $@".out" $^ $(OPTLIB)

clean:
	@rm *.o *.mod *.out
