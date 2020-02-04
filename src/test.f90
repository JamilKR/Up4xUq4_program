program test
  use MOD_Up4
  use MOD_matfun
  implicit none
  integer w,J,count
  !
  Npval=2
  count=0
  !
  do w=Npval,mod(Npval,2),-2
     !
     do J=0,w,1
        !
        write(*,10) "w1=",w,", J1=",J,", w2=",w,", J2=",J," ---> Qp=",RME_Qp2(w,J,w,J)
        write(*,10) "w1=",w,", J1=",J,", w2=",w+2,", J2=",J," ---> Qp=",RME_Qp2(w,J,w+2,J)
        write(*,10) "w1=",w,", J1=",J,", w2=",w-2,", J2=",J," ---> Qp=",RME_Qp2(w,J,w-2,J)
        !
        write(*,10) "w1=",w,", J1=",J,", w2=",w,", J2=",J+2," ---> Qp=",RME_Qp2(w,J,w,J+2)
        write(*,10) "w1=",w,", J1=",J,", w2=",w+2,", J2=",J+2," ---> Qp=",RME_Qp2(w,J,w+2,J+2)
        write(*,10) "w1=",w,", J1=",J,", w2=",w-2,", J2=",J+2," ---> Qp=",RME_Qp2(w,J,w-2,J+2)
        !
        write(*,10) "w1=",w,", J1=",J,", w2=",w,", J2=",J-2," ---> Qp=",RME_Qp2(w,J,w,J-2)
        write(*,10) "w1=",w,", J1=",J,", w2=",w+2,", J2=",J-2," ---> Qp=",RME_Qp2(w,J,w+2,J-2)
        write(*,10) "w1=",w,", J1=",J,", w2=",w-2,", J2=",J-2," ---> Qp=",RME_Qp2(w,J,w-2,J-2)

        count=count+1
        print*,count
        
     enddo
     !
  enddo

10 format(A,I2,A,I2,A,I2,A,I2,A,F15.4)

end program test
