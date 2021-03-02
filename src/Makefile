
SHELL = /bin/sh
BIN   = ../bin/


FC = gfortran
FFLAGS = -g -Wall -O3 -fcheck=all -cpp
#FFLAGS = -Wall -O3 -cpp -fopenmp

MODLIB = -I/usr/local/lib/lapack95/lapack95_modules
OPTLIB = -lblas -llapack /usr/local/lib/lapack95/lapack95.a

BuckPgm_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90  minuit-cern.f BuckProgram_fit.f90
BuckPgm_obj = $(addsuffix .o, $(basename ${BuckPgm_src}))

IntPgm_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90  minuit-cern.f IntProgram_fit.f90
IntPgm_obj = $(addsuffix .o, $(basename ${IntPgm_src}))

test_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 test.f90
test_obj = $(addsuffix .o, $(basename ${test_src}))

test_2_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 test_2.f90
test_2_obj = $(addsuffix .o, $(basename ${test_2_src}))

test_3_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90 test_3.f90
test_3_obj = $(addsuffix .o, $(basename ${test_3_src}))

test_4_src = MOD_matfun.f90   test_4.f90
test_4_obj = $(addsuffix .o, $(basename ${test_4_src}))

test_5_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90  test_5.f90
test_5_obj = $(addsuffix .o, $(basename ${test_5_src}))

test_6_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90  test_6.f90
test_6_obj = $(addsuffix .o, $(basename ${test_6_src}))


test_7_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90  test_7.f90
test_7_obj = $(addsuffix .o, $(basename ${test_7_src}))


test_8_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90  test_8.f90
test_8_obj = $(addsuffix .o, $(basename ${test_8_src}))


test_9_src = MOD_matfun.f90 MOD_Up4.f90 MOD_Uq4.f90 MOD_Up_x_Uq.f90  test_9.f90
test_9_obj = $(addsuffix .o, $(basename ${test_9_src}))

%.o: %.f90 
	$(FC) $(FFLAGS) -c $(MODLIB) -o $@ $<

%.o: %.f
	$(FC) -O3 -c  -o $@ $<

all: help BuckProgram_fit IntProgram_fit#test test_2 test_3 test_4 test_5 test_6

help:
	@echo
	@echo   "\e[42mDemasiadas expectativas...!\e[0m"
	@echo
	@echo   "\e[44mToo many expectations...!\e[0m"
	@echo

BuckProgram_fit: $(BuckPgm_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

IntProgram_fit: $(IntPgm_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)	


test: $(test_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_2: $(test_2_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_3: $(test_3_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_4: $(test_4_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_5: $(test_5_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_6: $(test_6_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_7: $(test_7_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_8: $(test_8_obj)
	$(FC) $(FFLAGS)  -o $(BIN)$@".out" $^ $(OPTLIB)

test_9: $(test_9_obj)
	$(FC) $(FFLAGS)  -o $@".out" $^ $(OPTLIB)

clean:
	@rm -f *.o *.mod *.out *~ $(BIN)"BuckProgram_fit.out"
