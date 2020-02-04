module MOD_matfun
  !
  implicit none
  !
  ! Usefull Math Functions
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
end module MOD_matfun
    
