program test_4
  use MOD_matfun
  implicit none
  !
  integer:: j1,j2,j3,l1,l2,l3
  double precision:: wg
  !
  write(*,*)wigner_6j(1,2,3,4,1,2)
  write(*,*)wigner_6j(1,2,3,4,2,3)
  write(*,*)wigner_6j(1,2,3,4,3,4)
  !
end program test_4
