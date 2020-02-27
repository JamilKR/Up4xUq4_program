program test_7
  use MOD_Up4,only: RME_np,Npval
  implicit none
  !
  integer:: w1, w2, j1, j2
  double precision:: np
  !
  Npval = 4
  !
  do w1=Npval,0,-1
     do w2=Npval,0,-1
        do j1=0,w1
           j2=j1
           np = RME_np(w1,j1,w2,j2,iprint = 2)
        enddo
     enddo
  enddo
  !
end program test_7
