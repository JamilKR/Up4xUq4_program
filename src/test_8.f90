program test_8
  use MOD_matfun
  implicit none
  !
  integer:: j1,l1,lam1, j2,l2,lam2
  double precision :: a

  j1=1
  l1=0
  lam1=1

  j2=2
  l2=0
  lam2=2
  
  a=wigner_9j(j1,l1,lam1, j2,l2,lam2, 0,1,1)
  print*, 'k=0, m=1, T=1', a
  
  a=wigner_9j(j1,l1,lam1, j2,l2,lam2, 2,1,1)
  print*, 'k=2, m=1, T=1', a
  
  a=wigner_9j(j1,l1,lam1, j2,l2,lam2, 2,0,2)
  print*, 'k=2, m=0, T=2', a
    
  a=wigner_9j(j1,l1,lam1, j2,l2,lam2, 2,2,2)
  print*, 'k=2, m=0, T=2', a
end program test_8
  
