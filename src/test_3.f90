program test_3
  use MOD_Up_x_Uq
  implicit none

  integer:: lambda_max, Nq, Np, i, total_para,total_ortho
  integer,allocatable:: dim_para(:),dim_ortho(:),ijk(:)
  integer,allocatable:: basis_para(:,:),basis_ortho(:,:)
  Np=34
  Nq=20

  
  
  lambda_max=5

  allocate(dim_para(0:lambda_max),dim_ortho(0:lambda_max),ijk(0:lambda_max))

  dim_para=0
  dim_ortho=0

  call dimension_po(Np,Nq,lambda_max,dim_para,dim_ortho)
  write(*,*) 'lambd: | ',(i,i=0,lambda_max)
  write(*,*) '---------------------------------------------------------------------------------'
  write(*,*) ' para: | ',(dim_para(i),i=0,lambda_max)
  write(*,*) 'ortho: | ',(dim_ortho(i),i=0,lambda_max)

  !

  total_para=sum(dim_para)
  total_ortho=sum(dim_ortho)

  allocate(basis_para(1:5,1:total_para),basis_ortho(1:5,1:total_ortho))


  call build_basis_po(Np,Nq,lambda_max,total_para,total_ortho,dim_para,dim_ortho,basis_para,basis_ortho)

  write(*,*)
  write(*,*)
  write(*,*) '---------------------------------------------------------------------------------'
  write(*,*) "S - w - J - n - L - lam"
  write(*,*) "-----------------------"
  do i=1,total_para
     write(*,10) "p -",basis_para(1,i),"-",basis_para(2,i),"-",basis_para(3,i),"-",basis_&
          &para(4,i),"-",basis_para(5,i)
  enddo
  write(*,*) "-----------------------"
  do i=1,total_ortho
     write(*,10) "o -",basis_ortho(1,i),"-",basis_ortho(2,i),"-",basis_ortho(3,i),"-",basi&
          &s_ortho(4,i),"-",basis_ortho(5,i)
  enddo

  
  write(*,*) '---------------------------------------------------------------------------------'

  call initialize_position_index(ijk,dim_ortho,lambda_max)
  print*,ijk(5),dim_ortho(5)
  do i = ijk(5),dim_ortho(5)+ijk(5)-1
     write(*,10) "p -",basis_ortho(1,i),"-",basis_ortho(2,i),"-",basis_ortho(3,i),"-",basis_&
          &ortho(4,i),"-",basis_ortho(5,i)
  enddo
  
10 format(5(A,I3))
end program test_3
