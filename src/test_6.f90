program test_6
  use MOD_matfun
  use MOD_Up4
  use MOD_Uq4
  use MOD_Up_x_Uq
  implicit none
  integer:: basis(1:5,1:10)
  double precision:: matrix(1:10,1:10)
  Npval=20
  Nqval=34

  basis=0
  
  call build_Up_x_Uq_matrix(basis=basis,matrix=matrix,RME_fun=RME_Ip_x_SOq4,iprint=3)
 
end program test_6
