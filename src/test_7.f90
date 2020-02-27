program test_7
  use MOD_Up4,only: RME_np,Npval
  use MOD_Uq4,only: RME_Dq_prima,Nqval
  implicit none
  !
  integer:: n1, n2, l1, l2
  double precision:: dq
  !
  Nqval = 4
  !
  do n1=Nqval,0,-1
     do n2=Nqval,0,-1
        do l1=n1,0,-2
           do l2=n2,0,-2
              dq = RME_Dq_prima(n1,l1,n2,l2,iprint = 2)
           enddo
        enddo
     enddo
  enddo
  !
end program test_7
