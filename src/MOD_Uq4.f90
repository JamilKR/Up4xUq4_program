module MOD_Uq4
  !
  implicit none
  !
  integer:: Nqval ! Uq(4) totally symmetric representation
  !
contains
  !
  function RME_Casimir_Uq3(n1,l1,n2,l2)
    !
    ! < [Nq] n1 l1 || nq || [Nq] n2 l2 >
    !
    implicit none
    !
    integer,intent(in):: n1,l1,n2,l2
    double precision:: RME_Casimir_Uq3
    !
    if ( (n1/=n2) .or. (l1/=l2) ) then
       !
       RME_Casimir_Uq3 = 0.0d0
       !
    else
       !
       RME_Casimir_Uq3 = dble(n1)
       !
    endif
    !
  end function RME_Casimir_Uq3
  !
  !*****************************************************************************************
  !
  function RME_Casimir_SOq3(n1,l1,n2,l2)
    !
    ! < [Nq] n1 l1 || L**2 || [Nq] n2 l2 >
    !
    implicit none
    !
    integer,intent(in):: n1,l1,n2,l2
    double precision:: RME_Casimir_SOq3
    !
    if ( (n1/=n2) .or. (l1/=l2) ) then
       !
       RME_Casimir_SOq3 = 0.0d0
       !
    else
       !
       RME_Casimir_SOq3 = dble(l1*(l1+1))
       !
    endif
    !
  end function RME_Casimir_SOq3
  !
  !*****************************************************************************************
  !
  function RME_Casimir_SOq4(n1,l1,n2,l2)
    !
    ! < [Nq] n1 L1 || C_2[SOq(4)] || [Nq] n2 L2 >
    ! < [Nq] n1 L1 || L**2 + D**2 || [Nq] n2 L2 >
    !
    implicit none
    !
    integer,intent(in):: n1,l1,n2,l2
    double precision:: RME_Casimir_SOq4
    !
    if (l1 /= l2) then
       !
       RME_Casimir_SOq4 = 0.0d0
       !
    else
       !
       if (n1 == n2) then
          !
          RME_Casimir_SOq4 = dble(Nqval*(2*n1+3)) - 2.0d0*dble(n1*(n1+1)) + dble(l1*(l1+1))
          !
       elseif (n1 == n2+2) then
          !
          RME_Casimir_SOq4 = sqrt(dble(Nqval-n2-1)*dble(Nqval-n2)*dble(n2-l2+2)*dble(n2+l2+3))
          !
       elseif (n1 == n2-2) then
          !
          RME_Casimir_SOq4 = sqrt(dble(Nqval-n2+1)*dble(Nqval-n2+2)*dble(n2-l2)*dble(n2+l2+1))
          !
       else
          !
          RME_Casimir_SOq4 = 0.0d0
          !
       endif
       !
    endif
    !
  end function RME_Casimir_SOq4
  !
  !*****************************************************************************************
  !
  function RME_Qq2(n1,l1,n2,l2)
    !
    ! < [Nq] n1 L1 || Qq || [Nq] n2 L2 >
    ! Corrected with factor to DS-Talmi.
    !
    implicit none
    !
    integer, intent(in):: n1,l1,n2,l2
    double precision:: RME_Qq2
    !
    if (n1 /= n2) then
       !
       RME_Qq2 = 0.0d0
       !
    else
       !
       if (l1 == l2) then
          !
          RME_Qq2 = dble(2*n2+3)*sqrt( (dble(2*l2+1)*dble(l2)*dble(l2+1)) / &
               (6.0d0*dble(2*l2-1)*dble(2*l2+3)) )
          !
       elseif (l1 == l2+2) then
          !
          RME_Qq2 = sqrt( (dble(n2-l2)*dble(n2+l2+3)*dble(l2+1)*dble(l2+2)) / dble(2*l2+3) )
          !
       elseif (l1 == l2-2) then
          !
          RME_Qq2 = sqrt( (dble(n2-l2+2)*dble(n2+l2+1)*dble(l2)*dble(l2-1)) / dble(2*l2-1) )
          !
       else
          !
          RME_Qq2 = 0.0d0
          !
       endif
       !
    endif
    !
  end function RME_Qq2
  !
  !*****************************************************************************************
  !
  function RME_Dq_prima(n1,l1,n2,l2,iprint)
    !
    ! < [Nq] n1 L1 || Dq' || [Nq] n2 L2 >
    !
    implicit none
    !
    integer, intent(in):: n1,l1,n2,l2
    double precision:: RME_Dq_prima
    integer, optional:: iprint
    !
    if (present(iprint) .and. iprint .ge. 2) then
       !
       write(*,111)  "< [",Nqval,"] ",n1,", ",l1," || Dq' || [",Nqval,"] ",n2,", ",l2," >"
111    format(A,2(I3,A,2(I2,A)))
       !
    endif
    !
    RME_Dq_prima = 0.0d0
    !
    if ( (n1 == n2-1) .and. (l1 == l2-1) ) then
       !
       RME_Dq_prima = - 1.0d0*sqrt( dble(Nqval-n2+1)*dble(n2+l2+1)*dble(l2) )
       !
    elseif ( (n1==n2-1) .and. (l1==l2+1) ) then
       !
       RME_Dq_prima = -1.0d0*sqrt( dble(Nqval-n2+1)*dble(n2-l2)*dble(l2+1) )
       !
    elseif ( (n1==n2+1) .and. (l1==l2-1) ) then
       !
       RME_Dq_prima = sqrt( dble(Nqval-n2)*dble(n2-l2+2)*dble(l2) )
       !
    elseif ( (n1==n2+1) .and. (l1==l2+1) ) then
       !
       RME_Dq_prima = sqrt( dble(Nqval-n2)*dble(n2+l2+3)*dble(l2+1) )
       !
    else
       !
       RME_Dq_prima = 0.0d0
       !
    endif
    !
    if(present(iprint) .and. iprint .ge. 2) write(*,"(T40,F7.2)") RME_Dq_prima
    !
  end function RME_Dq_prima
  !
  !*****************************************************************************************
  !
end module MOD_Uq4
