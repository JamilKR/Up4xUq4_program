module MOD_matfun
  !
  implicit none
  !
  ! Useful Math Functions
  !
contains
  !
  function p_symbol(a,s)
    !
    ! Pochhammer Symbol
    !
    ! (a)_s = a * (a+1) * ... * (a+s-1)
    !
    implicit none
    !
    integer,intent(in):: a,s
    double precision:: p_symbol
    integer:: i
    !
    p_symbol=dble(a)
    !
    do i=a+1,a+s-1
       !
       p_symbol = p_symbol*dble(i)
       !
    enddo
    !
  end function p_symbol
  !
  !*****************************************************************************************
  !
  function factorial(n)
    !
    ! No explanation is needed.
    !
    ! factorial(n) = n!
    !
    implicit none
    !
    integer, intent(in):: n
    double precision:: factorial
    integer:: i
    !
    if (n .lt. 0) STOP "ERROR! Are you trying to compute the factorial &
         &of a negative number?"
    !
    factorial=1.0d0
    !
    do i=2,n
       !
       factorial=factorial*dble(i)
       !
    end do
    !
  end function factorial
  !
  !*****************************************************************************************
  !
  function delta_function(a,b,c)
    !
    ! Function defined in eq. (8.63) of the book Lie Algebras and App. of A. Iachello.
    !
    implicit none
    !
    integer,intent(in):: a,b,c
    double precision:: delta_function
    !
    delta_function = sqrt( (factorial(a+b-c)*factorial(b+c-a)*factorial(c+a-b) ) / &
         factorial(a+b+c+1) )
    !
  end function delta_function
  !
  !*****************************************************************************************
  !
  function wigner_6j(j1,j2,j3,l1,l2,l3)
    !
    ! eq. (15.38) DS-Talmi
    ! Nuclear Shell Theory 
    !
    ! LaTeX: \{ j1 j2 j3 \\ l1 l2 l3 \}
    !
    implicit none
    !
    integer,intent(in):: j1,j2,j3,l1,l2,l3
    integer:: t_min, t_max, t
    double precision:: wigner_6j
    !
    wigner_6j = 0.0d0
    !
    ! Triangular Condition:
#define aaa j1
#define bbb j2
#define ccc j3
    if ( (abs(aaa+bbb) .lt. ccc) .or. (abs(aaa-bbb) .gt. ccc) ) then ! 1
       !
       return
       !
    end if
#undef aaa
#undef bbb
#undef ccc
    !
#define aaa j1
#define bbb l2
#define ccc l3
    if ( (abs(aaa+bbb) .lt. ccc) .or. (abs(aaa-bbb) .gt. ccc) ) then ! 2
       !
       return
       !
    end if
#undef aaa
#undef bbb
#undef ccc
    !
#define aaa l1
#define bbb j2
#define ccc l3
    if ( (abs(aaa+bbb) .lt. ccc) .or. (abs(aaa-bbb) .gt. ccc) ) then ! 3
       !
       return
       !
    end if
#undef aaa
#undef bbb
#undef ccc
    !
#define aaa l1
#define bbb l2
#define ccc j3
    if ( (abs(aaa+bbb) .lt. ccc) .or. (abs(aaa-bbb) .gt. ccc) ) then ! 4
       !
       return
       !
    end if
#undef aaa
#undef bbb
#undef ccc
    !
    t_min = max(j1+j2+j3,j1+l2+l3,l1+j2+l3,l1+l2+j3)
    t_max = min(j1+j2+l1+l2,j2+j3+l2+l3,j1+j3+l1+l3)
    !
    if (t_min .gt. t_max) STOP "ERROR! wigner_6j: t_min > t_max"
    !
    do t = t_min,t_max
       !
       wigner_6j = wigner_6j  + ( dble( (-1)**t ) )* ( factorial(t+1) / &
            ( factorial(t-j1-j2-j3) * factorial(t-j1-l2-l3) * factorial(t-l1-j2-l3) * &
            factorial(t-l1-l2-j3) * factorial(j1+j2+l1+l2-t) * factorial(j2+j3+l2+l3-t) * &
            factorial(j1+j3+l1+l3-t) ) ) 
       !
    enddo
    !
    wigner_6j = wigner_6j * delta_function(j1,j2,j3) * delta_function(j1,l2,l3) * &
         delta_function(l1,j2,l3)*delta_function(l1,l2,j3)
    !
  end function wigner_6j
  !
  !*****************************************************************************************
  !
end module MOD_matfun
    
