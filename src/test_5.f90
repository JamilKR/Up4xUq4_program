program test_5
  use MOD_matfun
  use MOD_Up4
  use MOD_Uq4
  use MOD_Up_x_Uq
  implicit none

  integer:: Np,Nq, total_para,total_ortho
  integer,parameter:: lambda_max=2
  integer:: dim_para(0:lambda_max),dim_ortho(0:lambda_max)!,ijk(0:lambda_max)
  integer,allocatable:: basis_para(:,:),basis_ortho(:,:)
  double precision, allocatable:: A(:,:), B(:,:)
  integer:: i,j
  Np=3
  Npval=3
  Nq=2
  Nqval=2

  dim_para=0
  dim_ortho=0

  call dimension_po(Np,Nq,lambda_max,dim_para,dim_ortho)
  
  total_para=sum(dim_para)
  total_ortho=sum(dim_ortho)
  allocate(basis_para(1:5,1:total_para),basis_ortho(1:5,1:total_ortho))

  
  call build_basis_po(Np,Nq,lambda_max,total_para,total_ortho, &
       dim_para,dim_ortho,basis_para,basis_ortho)
  
  allocate(A(1:total_para,1:total_para), B(1:total_para,1:total_para))

 ! A=build_Qp_x_Qq_0_matrix(basis_para,total_para)

  print*,"***********************************************************************"
  
  !B=build_Qp_x_Qq_0_matrix(basis_ortho,total_ortho)

  ! print*, total_para, total_ortho

!   write(*,*) "****** Para ******"
!   write(*,*) "******************"
!   write(*,*)

!   do i=1,total_para
!      write(*,10) (A(i,j) ,j=1,total_para)
!   enddo
! 10 format(31F6.1)
  
!   write(*,*) "****** Ortho ******"
!   write(*,*) "******************"
!   write(*,*)
  
!   do i=1,total_ortho
!      write(*,20) (B(i,j) ,j=1,total_ortho)
!   enddo
! 20 format(22F6.1)
  

end program test_5
