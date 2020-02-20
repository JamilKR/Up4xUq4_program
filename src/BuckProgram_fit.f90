program BuckProgram
  use MOD_matfun
  use MOD_Up4
  use MOD_Uq4
  use MOD_Up_x_Uq
#ifdef _OPENMP
  use omp_lib,only:omp_get_thread_num
#endif
  implicit none
  !
  character(len=50):: input, form, inp_exp !input namelist, format string, input exp trans
  integer:: i,j,k ! Loops
  integer:: ird=5, iwr =6, isav=7
  ! Hamiltonian parameters:
  double precision:: &
       ! H molecule:
       bet  , & ! C_2[SOp(4)]
       gam  , & ! C_2[SOp(3)]
       gam2 , & ! C_2[SOp(3)]**2
       kap  , & ! C_2[SOp(3)]*C_2[SOp(4)]
       ! H cage:
       a , & ! C_1[Uq(3)]
       b , & ! C_2[Uq(3)]
       c , & ! C_2[SOq(3)]
       d , & ! C_2[SOq(4)]
       ! H interaction:
       Qpq  , & ![Qp x Qq]^(0)
       Qpqw , & !C_2[SOp(4)][Qp x Qq]^(0) + [Qp x Qq]^(0)C_2[SOp(4)]
       v1       !C_1[Uq(3)]C_2[SOp(4)]
  !
  double precision,allocatable:: aux(:,:) ! Auxiliary matrix!!
  character(len=50) :: fit_par ! Parameters to be fixed
  !
  NAMELIST /num/ Npval, Nqval, lambda_max, inp_exp
  NAMELIST /mol/ bet, gam, gam2, kap
  NAMELIST /cag/ a, b, c, d
  NAMELIST /int/ Qpq, Qpqw, v1
  NAMELIST /fit/ Fit_par
  !
  write(*,'(A)') "Reading input file..."
  read(*,*) input
  open(unit=11,file=trim(input),status='old',action='read')
  !
  read(11,num)
  read(11,mol)
  read(11,cag)
  read(11,int)
  read(11,fit)
  !
  close(11)
  !
  allocate(dim_para(0:lambda_max), dim_ortho(0:lambda_max), &
       ijk_para(0:lambda_max),ijk_ortho(0:lambda_max))
  !
  write(*,'(/,A)') "Calculating dimensions..."
  call dimension_po(Npval,Nqval,lambda_max,dim_para,dim_ortho)
  !
  ! Format variable:
  write(form,*) "(T5,A,", lambda_max+1,"I10)"
  !
  write(*,form) "   Lambda =", (i, i=0,lambda_max)
  !
  write(*,form) " Para-dim =",(dim_para(i),i=0,lambda_max)
  write(*,form) "Ortho-dim =",(dim_ortho(i),i=0,lambda_max)
  !
  ! Allocate basis and compute them:
  allocate(basis_para(1:5,1:sum(dim_para)),basis_ortho(1:5,1:sum(dim_ortho)))
  call build_basis_po(Npval,Nqval,lambda_max,sum(dim_para),sum(dim_ortho), &
       dim_para,dim_ortho,basis_para,basis_ortho)
  !
  !
  write(*,'(/,A)') "Reading experimental energies and transitions..."
  total_exp=exp_lines(inp_exp,21)
  allocate(exp_data(1:total_exp))
  call read_expdat(inp_exp,exp_data,total_exp,21)
  !
  write(*,'(/,A)') "Building operator matrices..."
  !
  allocate(Ham(0:lambda_max),SOq4(0:lambda_max),QpQq(0:lambda_max),QpQqW(0:lambda_max))
  !
  call initialize_position_index(ijk_para,dim_para,lambda_max)
  call initialize_position_index(ijk_ortho,dim_ortho,lambda_max)
  !
  !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(aux) 
  do i=0,lambda_max
     !
     allocate( Ham(i)%para(1:dim_para(i),1:dim_para(i)), &
          SOq4(i)%para(1:dim_para(i),1:dim_para(i)), &
          QpQq(i)%para(1:dim_para(i),1:dim_para(i)), &
          QpQqW(i)%para(1:dim_para(i),1:dim_para(i)) )
     !
     allocate( Ham(i)%ortho(1:dim_ortho(i),1:dim_ortho(i)), &
          SOq4(i)%ortho(1:dim_ortho(i),1:dim_ortho(i)), &
          QpQq(i)%ortho(1:dim_ortho(i),1:dim_ortho(i)), &
          QpQqW(i)%ortho(1:dim_ortho(i),1:dim_ortho(i)) )
     !
     !ijk(i):ijk(i)+dim_para(i)-1 := from first state of the lambda-block to the last one
     ! SOq4 para-ortho
     call build_Up_x_Uq_matrix(basis_para(1:5,ijk_para(i):ijk_para(i)+dim_para(i)-1), &
          SOq4(i)%para, RME_Ip_x_SOq4)
     call build_Up_x_Uq_matrix(basis_ortho(1:5,ijk_ortho(i):ijk_ortho(i)+dim_ortho(i)-1), &
          SOq4(i)%ortho, RME_Ip_x_SOq4)
     !
     call build_Up_x_Uq_matrix(basis_para(1:5,ijk_para(i):ijk_para(i)+dim_para(i)-1), &
          QpQq(i)%para,RME_Qp_x_Qq_0)
     call build_Up_x_Uq_matrix(basis_ortho(1:5,ijk_ortho(i):ijk_ortho(i)+dim_ortho(i)-1), &
          QpQq(i)%ortho,RME_Qp_x_Qq_0)
     !
     Uncomment !!!! 
     allocate(aux(1:dim_para(i),1:dim_para(i)))
     call build_Up_x_Uq_matrix(basis_para(1:5,ijk_para(i):ijk_para(i)+dim_para(i)-1), &
          aux, RME_Ip_x_SOq4) ! SOp4 diag-matrix para base
     QpQqW(i)%para=0.5d0*( matmul(QpQq(i)%para,aux) + matmul(aux,QpQq(i)%para) )
     deallocate(aux)
     !
     allocate(aux(1:dim_ortho(i),1:dim_ortho(i)))
     call build_Up_x_Uq_matrix(basis_ortho(1:5,ijk_ortho(i):ijk_ortho(i)+dim_ortho(i)-1),&
          aux, RME_Ip_x_SOq4) ! SOp4 diag-matrix ortho base
     QpQqW(i)%ortho=0.5d0*( matmul(QpQq(i)%ortho,aux) + matmul(aux,QpQq(i)%ortho) )
     deallocate(aux)
     !
#ifdef _OPENMP
     write(*,'(T5,A,I3,A,I3)') "Computed matrices for lambda =", i," on thread ", &
          omp_get_thread_num()
#else
     write(*,'(T5,A,I3)') "Computed matrices for lambda =", i
#endif
  enddo
  !$OMP END PARALLEL DO
  !
  open(unit=12,file='minuit_input.inp',status='replace')
  inquire(12)
  !
  write(12,'(A)') "SET TITLE"
  write(12,'(A)') "'Minuit minimization'"
  write(12,'(A)') "PARAMETERS"
  write(12,'(A,D20.10,A)') "1    'bet '",  bet,   "         0.1D-02"
  write(12,'(A,D20.10,A)') "2    'gam '",  gam,   "         0.1D-02"
  write(12,'(A,D20.10,A)') "3    'gam2 '", gam2,  "         0.1D-02"
  write(12,'(A,D20.10,A)') "4    'kap '",  kap,   "         0.1D-02"
  write(12,'(A,D20.10,A)') "5    'a '",    a,     "         0.1D-02"
  write(12,'(A,D20.10,A)') "6    'b '",    b,     "         0.1D-02"
  write(12,'(A,D20.10,A)') "7    'c '",    c,     "         0.1D-02"
  write(12,'(A,D20.10,A)') "8    'd '",    d,     "         0.1D-02"
  write(12,'(A,D20.10,A)') "9    'Qpq '",  Qpq,   "         0.1D-02"
  write(12,'(A,D20.10,A)') "10   'Qpqw '", Qpqw,  "         0.1D-02"
  write(12,'(A,D20.10,A)') "11   'v1 '",   v1,    "         0.1D-02"
  !
  write(12,*)
  write(12,'(A)') adjustl(fit_par)
  write(12,'(A)') "set stra 2"
  write(12,'(A)') "minimize 10000"
  write(12,'(A)') "call 3"
  write(12,'(A)') "exit"
  write(12,*)
  write(12,*)
  !
  close(12)
  !
  write(*,'(/,A)') "Calling minuit ..."
  !
  !params = (/ bet, gam, gam2, kap, a, b, c, d, Qpq, Qpqw, v1 /)
  open(unit=ird,file='minuit_input.inp',status='old')
  call mintio(ird,iwr,isav) ! units read, write and save
  call minuit(FCN,chi2)
  !
  !
end program BuckProgram
