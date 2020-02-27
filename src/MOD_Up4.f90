module MOD_Up4
  !
  use MOD_matfun,only: p_symbol
  !
  implicit none
  !
  integer:: Npval ! Up(4) totally symmetric representation
  !
contains
  !
  function RME_Casimir_SOp4(w1,j1,w2,j2)
    !
    ! < [Np] w1 j1 || C2[sop(4)] || [Np] w2 j2 >
    !
    implicit none
    !
    integer,intent(in):: w1,j1,w2,j2
    double precision:: RME_Casimir_SOp4
    integer:: v
    !
    if ((w1/=w2) .or. (j1/=j2)) then
       !
       RME_Casimir_SOp4 = 0.0d0
       !
    else
       !
       !RME_Casimir_SOp4 = dble( w1*(w1+2) )
       v = (Npval-w1)/2
       RME_Casimir_SOp4 = 4.0d0*dble(v)*(1.0d0/dble(Npval)-1.0d0+dble(v)/dble(Npval))
       !
    endif
    !
  end function RME_Casimir_SOp4
  !
  !*****************************************************************************************
  !
  function RME_Casimir_SOp3(w1,j1,w2,j2)
    !
    ! < [Np] w1 j1 || C2[sop(4)] || [Np] w2 j2 >
    !
    implicit none
    !
    integer,intent(in):: w1,j1,w2,j2
    double precision:: RME_Casimir_SOp3
    !
    if ((w1/=w2) .or. (j1/=j2)) then
       !
       RME_Casimir_SOp3 = 0.0d0
       !
    else
       !
       RME_Casimir_SOp3 = dble( j1*(j1+1) )
       !
    endif
    !
  end function RME_Casimir_SOp3
  !
  !*****************************************************************************************
  !
  function RME_Qp2(w1,j1,w2,j2,iprint)
    !
    ! < [Np] w1 j1 || Qp2 || [Np] w2 j2 >
    !
    ! Corrected with factor to DS-Talmi.
    !
    implicit none
    !
    integer, intent(in):: w1,j1,w2,j2
    double precision:: RME_Qp2
    integer, optional:: iprint
    integer:: iprint2
    !
    if (present(iprint)) then
       iprint2 = iprint
    else
       iprint2 = 0
    endif
    !
    if (j1 == j2) then
       !
       if(iprint2 .ge. 1) write(*,'(A,I3)') "RME_Qp2 j1 = j2   = ", j1
       !
       if (w1 == w2) then ! < w J || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2   = ", w1
          !
          if (j2 == 0) then
             !
             RME_Qp2 = 0.0d0
             !
          else
             !
             RME_Qp2 = dble(Npval+2) * (1.0d0+dble(j2*(j2+1))/dble(w2*(w2+2))) * &
                  sqrt( dble(j2*(j2+1))*dble((2*j2+1))/ &
                  ( 6.0d0*dble((2*j2-1))*dble((2*j2+3)) ) )
             !
          endif
          !
       elseif (w1 == w2+2) then ! < w+2 J || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2+2 = ", w1
          !
          RME_Qp2 = sqrt( ( dble(Npval-w2)*dble(Npval+w2+4)*p_symbol(w2-j2+1,2)* &
               p_symbol(w2+j2+2,2)*dble(j2*(j2+1))*dble((2*j2+1)) ) / &
               ( 24.0d0*p_symbol(w2+1,3)*dble(w2+2)*dble(2*j2-1)*dble(2*j2+3) )  )
          !
       elseif (w1 == w2-2) then ! < w-2 J || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2-2 = ", w1
          !
          RME_Qp2 = sqrt( ( dble(Npval-w2+2)*dble(Npval+w2+2)*p_symbol(w2-j2-1,2)* &
               p_symbol(w2+j2,2)*dble(j2*(j2+1))*dble((2*j2+1)) ) / &
               ( 24.0d0*p_symbol(w2-1,3)*dble(w2)*dble(2*j2+3)*dble(2*j2-1) ) )
          !
       else
          !
          if(iprint2 .ge. 1) write(*,'(A)') "        w1 = ELSE"
          !
          RME_Qp2 = 0.0d0
          !
       endif
       !
    elseif (j1 == j2+2) then
       !
       if(iprint2 .ge. 1) write(*,'(A,I3)') "RME_Qp2 j1 = j2+2 = ", j1
       !
       if (w1 == w2) then ! < w J+2 || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2   = ", w1
          !
          RME_Qp2 = dble(Npval+2)*sqrt( ( p_symbol(w2-j2-1,2)*p_symbol(w2+j2+2,2)* &
               dble(j2+2)*dble(j2+1) ) / ( 4.0d0*dble(w2*(w2+2))*dble(w2*(w2+2))* &
               dble(2*j2+3) ) )
          !
       elseif (w1 == w2+2) then ! < w+2 J+2 || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2+2 = ", w1
          !
          RME_Qp2 = sqrt( ( dble(Npval-w2)*dble(Npval+w2+4)*p_symbol(w2+j2+2,4)* &
               dble(j2+2)*dble(j2+1) ) / ( 16.0d0*p_symbol(w2+1,3)*dble(w2+2)* &
               dble(2*j2+3) ) )
          !
       elseif (w1 == w2-2) then ! < w-2 J+2 || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2-2 = ", w1
          !
          RME_Qp2 = sqrt( ( dble(Npval-w2+2)*dble(Npval+w2+2)*p_symbol(w2-j2-3,4)* &
               dble(j2+2)*dble(j2+1) ) / ( 16.0d0*p_symbol(w2-1,3)*dble(w2)*dble(2*j2+3) ) )
          !
       else
          !
          if(iprint2 .ge. 1) write(*,'(A)') "        w1 = ELSE"
          !
          RME_Qp2 = 0.0d0
          !
       endif
       !
    elseif (j1 == j2-2) then
       !
       if(iprint2 .ge. 1) write(*,'(A,I3)') "RME_Qp2 j1 = j2-2 = ", j1
       !
       if (w1 == w2) then ! < w J-2 || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2   = ", w1
          !
          RME_Qp2 = dble(Npval+2)*sqrt( ( p_symbol(w2-j2+1,2)*p_symbol(w2+j2,2)* &
               dble(j2)*dble(j2-1) ) / ( 4.0d0*dble(w2*(w2+2))*dble(w2*(w2+2))* &
               dble(2*j2-1) ) )
          !
       elseif (w1 == w2+2) then ! < w+2 J-2 || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2+2 = ", w1
          !
          RME_Qp2 = sqrt( ( dble(Npval-w2)*dble(Npval+w2+4)*p_symbol(w2-j2+1,4)* &
               dble(j2)*dble(j2-1) ) / ( 16.0d0*p_symbol(w2+1,3)*dble(w2+2)*dble(2*j2-1) ) )
          !
       elseif (w1 == w2-2) then ! < w-1 J-2 || Qp || w J >
          !
          if(iprint2 .ge. 1) write(*,'(A,I3)') "        w1 = w2-2 = ", w1
          !
          RME_Qp2 = sqrt( ( dble(Npval-w2+2)*dble(Npval+w2+2)*p_symbol(w2+j2-2,4)* &
               dble(j2)*dble(j2-1) ) / ( 16.0d0*p_symbol(w2-1,3)*dble(w2)*dble(2*j2-1) ) )
          !
       else
          !
          if(iprint2 .ge. 1) write(*,'(A)') "        w1 = ELSE"
          !
          RME_Qp2 = 0.0d0
          !
       endif
       !
    else
       !
       if(iprint2 .ge. 1) write(*,'(A)') "RME_Qp2 j1 = ELSE "
       !
       RME_Qp2 = 0.0d0
       !
    endif
    !
  end function RME_Qp2
  !
  !*****************************************************************************************
  !
  function RME_np(w1,j1,w2,j2,iprint)
    !
    ! < [Np] w1 J1 || np || [Np] w2 J2 >
    !
    implicit none
    !
    integer, intent(in):: w1,j1,w2,j2
    double precision:: RME_np
    integer, optional:: iprint
    !
    if (present(iprint) .and. iprint .ge. 2) then
       !
       write(*,111)  "< [",Npval,"] ",w1,", ",J1," || np || [",Npval,"] ",w2,", ",J2," >"
111    format(A,2(I3,A,2(I2,A)))
       !
    endif
    !
    RME_np = 0.0d0
    !
    if ( j1/=j2 ) then
       !
       if(present(iprint) .and. iprint .ge. 2) write(*,"(T40,F7.2)") RME_np
       return
       !
    endif
    !
    if ( w1 == w2 ) then
       !
       if ( w1 == 0 ) then
          RME_np = 3.0d0*dble(Npval)/4.0d0
       else
          RME_np = ( 0.5d0*dble(Npval-1) + dble(Npval+2)*0.5d0*dble(j1*(j1+1))/ &
               dble(w2*(w2+2)) ) * sqrt(2.0d0*dble(j1)+1.0d0)
       endif
       !
    elseif ( w1 == w2+2 ) then
       !
       RME_np = sqrt( dble(2*j1+1)*dble(Npval-w2)*dble(Npval+w2+4)*dble(w2-j1+1)* &
            dble(w2-j1+2)*dble(w2+j1+2)*dble(w2+j1+3) / ( 16.0d0*dble((w2+2)*(w2+1))* &
            dble((w2+2)*(w2+3)) ) )
       !
    elseif ( w1 == w2-2 ) then
       !
       RME_np = sqrt( dble(2*j1+1)*dble(Npval-w2+2)*dble(Npval+w2+2)*dble(w2-j1-1)* &
            dble(w2-j1)*dble(w2+j1)*dble(w2+j2+1) / ( 16.0d0*dble(w2*(w2-1))* &
            dble(w2*(w2+1)) ) )
       !
    endif
    !
    if(present(iprint) .and. iprint .ge. 2) write(*,"(T40,F7.2)") RME_np
    !
  end function RME_np
  !
  !*****************************************************************************************
  !
end module MOD_Up4
