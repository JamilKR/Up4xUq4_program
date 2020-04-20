program test
  use MOD_Uq4
  use MOD_matfun
  use f95_lapack,only: la_syevr
  implicit none
  integer::n1,n2,k1,k2,stat
  double precision,allocatable:: work(:)
  integer,parameter:: L= 2
  integer:: dim
  double precision,allocatable:: A(:,:),B(:)
  !
  Nqval=10
  dim=(Nqval-L -mod(L,2))/2 +1
  !
  allocate(A(1:dim,1:dim),B(1:dim),work(6*dim-1) )
  A=0.0d0
  B=0.0d0
  !
  k1=1
  !
  do n1=L,Nqval,2
     !
     k2=1
     !
     do n2=L,Nqval,2
        !
        A(k1,k2) = RME_Casimir_SOq4(n1,L,n2,L)
        !
        k2=k2+1
        !
     enddo
     !
     write(*,'(6F10.4)') (A(k1,k2),k2=1,dim)
     !
     k1=k1+1
     !
  enddo
  !
  ! Diagonalization:
  !
  ! call dsyev('N','U',dim,A,dim+1,B,work,3*dim-1,info)
  ! if (info/=0)then
  !    print*, "info =", info
  !    stop "todo mal"
  ! endif
  call la_syevr(A=A,W=B,JOBZ='V',UPLO='L')
  !
  write(*,*) "L=",L,": ", (B(k1),k1=1,dim) 
  !
  !
  if (allocated(work)) deallocate(work,STAT=stat)
  print*,'Work',stat
  ! if (allocated(A)) deallocate(A,STAT=stat)
  ! print*, 'A', stat
  !
  !
end program test


