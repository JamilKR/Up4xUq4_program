program test_8
  use MOD_matfun
  implicit none
  !
  double precision :: a
  
  a=wigner_9j(0,1,1, 0,1,1, 0,1,1)
  print*, a
  a=wigner_9j(0,0,0, 0,2,2, 0,1,1)
  print*, a
  a=wigner_9j(1,1,0, 1,1,1, 0,1,1)
  print*, a

  a=wigner_9j(0,1,1, 0,1,1, 2,1,1)
  print*, a
  a=wigner_9j(0,0,0, 0,2,2, 2,1,1)
  print*, a
  a=wigner_9j(1,1,0, 1,1,1, 2,1,1)
  print*, a

  
end program test_8
  
