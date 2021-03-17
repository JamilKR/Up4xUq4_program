program BuckSystem
  use MOD_matfun
  use MOD_Up4
  use MOD_Uq4
  use MOD_Up_x_Uq
#ifdef _OPENMP
  use omp_lib,only:omp_get_thread_num
#endif
  implicit none
  !
  character(len=50):: input, form !input namelist, format string
  integer:: i! Loops
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
  NAMELIST /num/ Npval, Nqval, lambda_max
  NAMELIST /mol/ bet, gam, gam2, kap
  NAMELIST /cag/ a, b, c, d
  NAMELIST /int/ Qpq, Qpqw, v1
  !
  write(*,'(A)') "Reading input file..."
  read(*,*) input
  open(unit=11,file=trim(input),status='old',action='read')
  !
  read(11,num)
  read(11,mol)
  read(11,cag)
  read(11,int)
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
  allocate(basis_para(1:5,1:sum(dim_para)), &
       basis_ortho(1:5,1:sum(dim_ortho)), &
       EnergiesPara(1:sum(dim_para)), &
       EnergiesOrtho(1:sum(dim_ortho)) )
  !
  call build_basis_po(Npval,Nqval,lambda_max,sum(dim_para),sum(dim_ortho), &
       dim_para,dim_ortho,basis_para,basis_ortho)
  !
  write(*,'(/,A)') "Building operator matrices..."
  !
  allocate(Ham(0:lambda_max),SOq4(0:lambda_max),QpQq(0:lambda_max),QpQqW(0:lambda_max))
  !
  call initialize_position_index(ijk_para,dim_para,lambda_max)
  call initialize_position_index(ijk_ortho,dim_ortho,lambda_max)
  !
  !$OMP PARALLEL DO DEFAULT(SHARED)
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
     call build_Up_x_Uq_matrix(basis_para(1:5,ijk_para(i):ijk_para(i)+dim_para(i)-1), &
          QpQqW(i)%para,RME_QpQqSOp4)
     call build_Up_x_Uq_matrix(basis_ortho(1:5,ijk_ortho(i):ijk_ortho(i)+dim_ortho(i)-1), &
          QpQqW(i)%ortho,RME_QpQqSOp4)
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
  ! Compute EigenSystem
  call eigensystem(EnergiesPara, EnergiesOrtho, lambda_max, &
       bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1, &
       .false.) ! .false. = No eigenvector output
  !
  ! Save Energies
  call save_energies(EnergiesPara, EnergiesOrtho, lambda_max )
  !
end program BuckSystem
  
  
