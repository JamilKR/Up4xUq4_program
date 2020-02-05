program test_2
  !
  use MOD_Uq4
  !
  implicit none
  !
  double precision, allocatable:: A(:,:)
  integer,parameter:: n=10
  integer:: l1,l2,dim, k1, k2
  !
  Nqval = 10
  dim = (n - mod(n,2))/2+1
  !
  allocate(A(1:dim,1:dim))
  !
  k1=1
  do l1=n,mod(n,2),-2
     k2=1
     do l2=n,mod(n,2),-2
        
        !
        A(k1,k2)=RME_Qq2(n,l1,n,l2)
        !
        k2=k2+1
     enddo
     !
     write(*,'(6F10.4)') (A(k1,k2),k2=1,dim)
     !
     k1=k1+1
     !
  enddo
  !
end program test_2

