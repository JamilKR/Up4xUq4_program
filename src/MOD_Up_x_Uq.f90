module MOD_Up_x_Uq
  !
  use MOD_matfun
  use MOD_Up4
  use MOD_Uq4
#ifdef __GFORTRAN__
  !Lapack95
  USE F95_LAPACK, ONLY: LA_SYEVR
#elif defined __INTEL_COMPILER
  USE LAPACK95,ONLY: SYEVR
#endif
  !
  implicit none
  !
  integer, allocatable:: basis_para(:,:),basis_ortho(:,:), dim_para(:), dim_ortho(:), &
       ijk_para(:), ijk_ortho(:)
  integer::lambda_max
  type exp_point
     !
     integer:: ist(1:5) ! initial state
     integer:: i_pos    ! position in para/ortho basis
     integer:: fst(1:5) ! final state
     integer:: f_pos    ! fosition in para/ortho basis
     double precision:: energy
     double precision:: intensity
     !
  end type exp_point
  integer:: total_exp ! number of experimental data
  type(exp_point),allocatable:: exp_data(:) ! experimental data structures
  !
  type matrix
     !
     double precision,allocatable:: para(:,:) ! para-matrix
     double precision,allocatable::ortho(:,:) ! ortho-matrix
     !
  end type matrix
  ! matrices:
  type(matrix), allocatable:: &
       Ham(:)    , & ! Hamiltonian para and ortho
       SOq4(:)   , & ! SOq(4) operator matrices para and ortho
       QpQq(:)   , & ! Qp x Qq operator matrices para and ortho
       QpQqW(:)      ! Qp x Qq SOp(4) matrices para and ortho

  double precision, allocatable:: &
       EnergiesPara(:), &
       EnergiesOrtho(:) ! Fixed energies to compute intensities
  !
  ! Array with the expected value of transitions:
  ! shape = len(exp data)x3
  double precision, allocatable:: intop(:,:)
  double precision:: temp ! Reduced temperature (kB*T)
  !
contains
  !
  subroutine dimension_po(Np,Nq,lambda_max,dim_para,dim_ortho)
    !
    ! INPUTs:
    !        o) Np: Up(4)
    !        o) Nq: Uq(4)
    !        o) lambda_max: maximum lambda to be considered
    !
    ! OUTPUTs:
    !         o) dim_para (0:lambda_max): dimensions of the diferent  para-blocks
    !         o) dim_ortho(0:lambda_max): dimensions of the diferent ortho-blocks
    !
    implicit none
    !
    integer,intent(in):: Np,Nq,lambda_max                               ! inputs
    integer,intent(out)::dim_para(0:lambda_max),dim_ortho(0:lambda_max) ! outputs
    integer::w,j,n,l,lam !loop control
    !
    ! para-case:
    !
    dim_para=0
    !
    J_p:do j = 0,Np,2 ! para J angular momentum
       L_p:do l = 0,Nq,1
          !
          CONDITION_p:if (abs(j-l) .le. lambda_max) then
             !
             LAMB_p: do lam = abs(j-l), min(j+l,lambda_max), 1
                !
                do w = Np,j,-2 ! para-case --> J even
                   do n = L,Nq,2
                      !
                      dim_para(lam) = dim_para(lam) + 1
                      !
                   enddo
                enddo
                !
             enddo LAMB_P
             !
          endif CONDITION_P
          !
       enddo L_P
    enddo J_P
    !
    ! ********************
    !
    ! ortho-case:
    !
    dim_ortho=0
    !
    J_o:do j = 1,Np,2 ! ortho J angular momentum
       L_o:do l = 0,Nq,1
          !
          CONDITION_o:if (abs(j-l) .le. lambda_max) then
             !
             LAMB_o: do lam = abs(j-l), min(j+l,lambda_max), 1
                !
                do w = Np,j,-2 ! ortho-case --> J odd
                   do n = L,Nq,2
                      !
                      dim_ortho(lam) = dim_ortho(lam) + 1
                      !
                   enddo
                enddo
                !
             enddo LAMB_O
             !
          endif CONDITION_O
          !
       enddo L_O
    enddo J_O
    !
  end subroutine dimension_po
  !
  !*****************************************************************************************
  !
  subroutine build_basis_po(Np,Nq,lambda_max,total_para,total_ortho, &
       dim_para,dim_ortho,basis_para,basis_ortho)
    !
    ! INPUTs:
    !        o) Np: Up(4)
    !        o) Nq: Uq(4)
    !        o) lambda_max: maximum lambda to be considered
    !        o) total_para:  total dimension of para  case
    !        o) total_ortho: total dimension of ortho case
    !        o) dim_para:  partial dimensions of  para-blocks
    !        o) dim_ortho: partial dimensions of ortho-blocks
    !
    ! OUTPUTs:
    !         o) basis_para (1:5,1:total_para): basis  para-blocks
    !         o) basis_ortho(1:5,1:total_ortho): basis ortho-blocks
    !
    implicit none
    !
    integer,intent(in):: Np,Nq,lambda_max,total_para,total_ortho  ! inputs
    integer,intent(in):: dim_para(0:lambda_max),dim_ortho(0:lambda_max)
    integer,intent(out):: basis_para(1:5,1:total_para), &
         basis_ortho(1:5,1:total_ortho)                           !outputs
    integer:: w,j,n,l,lam !loop control
    integer:: ijk(0:lambda_max)
    !
    ! para-case:
    !
    basis_para=0
    !
    ! Each lambda block starts in the next position after the end of the block before
    call initialize_position_index(ijk,dim_para,lambda_max)
    !
    !
    J_p:do j = 0,Np,2 ! para J angular momentum
       L_p:do l = 0,Nq,1
          !
          CONDITION_p:if (abs(j-l) .le. lambda_max) then
             !
             LAMB_p: do lam = abs(j-l), min(j+l,lambda_max), 1
                !
                do w = Np,j,-2 ! para-case --> J even
                   do n = L,Nq,2
                      !
                      basis_para(1,ijk(lam)) = w
                      basis_para(2,ijk(lam)) = J
                      basis_para(3,ijk(lam)) = n
                      basis_para(4,ijk(lam)) = L
                      basis_para(5,ijk(lam)) = lam
                      !
                      ijk(lam) = ijk(lam)+1
                      !
                   enddo
                enddo
                !
             enddo LAMB_P
             !
          endif CONDITION_P
          !
       enddo L_P
    enddo J_P
    !
    ! ********************
    !
    ! ortho-case:
    !
    basis_ortho=0
    !
    ! Each lambda block starts in the next position after the end of the block before
    call initialize_position_index(ijk,dim_ortho,lambda_max)
    !
    !
    J_o:do j = 1,Np,2 ! ortho J angular momentum
       L_o:do l = 0,Nq,1
          !
          CONDITION_o:if (abs(j-l) .le. lambda_max) then
             !
             LAMB_o: do lam = abs(j-l), min(j+l,lambda_max), 1
                !
                do w = Np,j,-2 ! ortho-case --> J odd
                   do n = L,Nq,2
                      !
                      basis_ortho(1,ijk(lam)) = w
                      basis_ortho(2,ijk(lam)) = J
                      basis_ortho(3,ijk(lam)) = n
                      basis_ortho(4,ijk(lam)) = L
                      basis_ortho(5,ijk(lam)) = lam
                      !
                      ijk(lam) = ijk(lam)+1
                      !
                   enddo
                enddo
                !
             enddo LAMB_O
             !
          endif CONDITION_O
          !
       enddo L_O
    enddo J_O
    !
  end subroutine build_basis_po
  !
  !*****************************************************************************************
  !
  subroutine initialize_position_index(ijk,partial_dim,lambda_max)
    !
    ! Inputs:
    !        o) partial_dim: Integer array dimension 0:lambda_max
    !        o) lambda_max
    !
    ! Output:
    !        o) ijk: Integer array dimension 0:lambda_max
    !
    ! This functions initializes the initial integer "pointer" ijk(0:lambda_max),
    ! so that ijk(0) = 1
    !         ijk(1) = position where lambda=1 block starts
    !         ...
    !         ijk(lambda_max) = position where lambda=lambda_max block stats
    !
    ! This function is going to be of vital importance during the program's development
    !
    implicit none
    !
    integer, intent(in) :: lambda_max,partial_dim(0:lambda_max)
    integer, intent(out):: ijk(0:lambda_max)
    integer:: i
    !
    ijk=1 ! All position are initialized, but only ijk(0) is ready
    !
    do i=1,lambda_max
       ijk(i)=ijk(i-1)+partial_dim(i-1)
    enddo
    !
  end subroutine initialize_position_index
  !
  !*****************************************************************************************
  !
  function RME_Qp_x_Qq_0(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    ! < w1 J1, n1 L1, lam1 || [ Qp^2 x Qq^2 ]^0 || w2 J2, n2 L2, lam2 >
    !
    ! sqrt(5) factor omited!
    !
    implicit none
    !
    integer,intent(in):: w1,j1,n1,l1,lam1, w2,j2,n2,l2,lam2
    double precision  :: RME_Qp_x_Qq_0
    !
    RME_Qp_x_Qq_0 = 0.0d0
    !
    if ( lam1 /= lam2 ) return
          !
    RME_Qp_x_Qq_0 = dble( (-1)**(l1+lam1+j2) ) * &
         wigner_6j(j1,l1,lam1, l2,j2,2) * &
         RME_Qq2(n1,l1,n2,l2) * &
         RME_Qp2(w1,j1,w2,j2)
    !
  end function RME_Qp_x_Qq_0
  !
  !*****************************************************************************************
  !
  function RME_Ip_x_SOq4(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    ! < w1 J1, n1 L1, lam1 || [ Ip x C_2[SOq(4)] ]^0 || w2 J2, n2 L2, lam2 >
    !
    implicit none
    !
    integer, intent(in):: w1,J1,n1,L1,lam1, w2,J2,n2,L2,lam2
    double precision:: RME_Ip_x_SOq4
    !
    RME_Ip_x_SOq4 = 0.0d0
    !
    if ( (lam1/=lam2) .or. (w1/=w2) .or. (j1/=j2) ) return ! Kronecker deltas
    !
    RME_Ip_x_SOq4 = dble( (-1)**(l1+lam1+j1) ) * sqrt( dble(2*lam1+1) ) * &
         wigner_6j(j1,l1,lam1,l2,j1,0) * RME_Casimir_SOq4(n1,l1,n2,l2)
    !
  end function RME_Ip_x_SOq4
  !
  !*****************************************************************************************
  !
  function RME_QpQqSOp4(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    ! It's not the most smart way, but It's useful
    ! < state 1 | | 0.5 * ( SOp(4)[QpxQq] + [QpxQq]SOp(4) ) | | state 2 > 
    !
    implicit none
    !
    integer, intent(in):: w1,J1,n1,L1,lam1, w2,J2,n2,L2,lam2
    double precision:: RME_QpQqSOp4
    !
    RME_QpQqSOp4 = 0.0d0
    !
    if ( lam1 /= lam2 ) return
    !
    RME_QpQqSOp4 = 0.5d0 * ( RME_Casimir_SOp4(w1,j1,w1,j1) + &
         RME_Casimir_SOp4(w2,j2,w2,j2) ) * RME_Qp_x_Qq_0(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
  end function RME_QpQqSOp4
  !
  !*****************************************************************************************
  !
  subroutine build_Up_x_Uq_matrix(basis,matrix,RME_fun,iprint)
    !
    ! This function build the para / ortho matrices using the given basis.
    ! This procedure can be used to build Up4 x Uq4 operators' matrices.
    !
    ! INPUTs:
    !        o) basis:   para/ortho basis
    !        o) matrix:  square matriz len(basis) x len(basis)
    !        o) RME_fun: function with w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2 dependences.
    !        o) iprint:  printing control
    !
    ! OUTPUT:
    !        o) matrix
    !
    ! All position corresponding to lam1 /= lam2 will be ZERO! 
    !
    implicit none
    !
    integer,intent(in)           :: basis(:,:)
    double precision,intent(out) :: matrix(:,:)
    double precision,external    :: RME_fun
    integer,optional             :: iprint
    integer:: &
         lbdim(1:2), & ! Basis lower limit
         ubdim(1:2), & ! Basis upper limit
         lmtx(1:2),  & ! Matrix lower limits
         umtx(1:2),  & ! Matrix upper limits
         iprint2,    & ! If iprint is not given
         i,j
    !
    if ( present(iprint) ) then 
       iprint2 = iprint
    else
       iprint2 = 0
    endif
    !
    lbdim = lbound(basis)
    ubdim = ubound(basis)
    !
    lmtx = lbound(matrix)
    umtx = ubound(matrix)
    !
    ! Checking dimensions:
    !
    if ( (lbdim(1)/=1) .or. (ubdim(1)/=5) ) STOP "ERROR! build_Up_x_Uq_matrix: &
         &5 quantum numbers are needed !!! "
    !
    if ( (lmtx(1)/=lmtx(2)) .or. (umtx(1)/=umtx(2)) ) STOP "ERROR! build_Up_x_Uq_matrix: &
         &matrix must be a square matrix !!! "
    !
    if ( (umtx(1)/=ubdim(2)) .or. (lmtx(1)/=lbdim(2)) ) STOP "ERROR! build_Up_x_Uq_matrix: &
         &the dimensions of the basis and the matrix are different !!! "
    !
    matrix = 0.0d0
    !
    do i = lmtx(1),umtx(1) ! bra-index
       !
       if (iprint2 .ge. 1) write(*,"(/,A,/)") trim(pretty_braket(w=basis(1,i), &
            j=basis(2,i),n=basis(3,i),l=basis(4,i),lam=basis(5,i),bk='b',Np=Npval,Nq=Nqval))
       !
       do j = lmtx(2),umtx(2) ! ket-index
          !
          if (iprint2 .ge. 2) write(*,"(T30,A)") trim(pretty_braket(w=basis(1,j), &
               j=basis(2,j),n=basis(3,j),l=basis(4,j),lam=basis(5,j),bk='k', &
               Np=Npval,Nq=Nqval))
          !
          if ( basis(5,i) /= basis(5,j) ) cycle ! lambda_i must be .eq. to lambda_j
          !
          matrix(i,j) = RME_fun(basis(1,i),basis(2,i),basis(3,i),basis(4,i),basis(5,i), &
               basis(1,j),basis(2,j),basis(3,j),basis(4,j),basis(5,j))
          !
          if (iprint2 .ge. 3) write(*,'(/,T50,ES40.5E5,/)') matrix(i,j)
          !
       enddo
       !
    enddo
    !
  end subroutine build_Up_x_Uq_matrix
  !
  !*****************************************************************************************
  !
  function pretty_braket(w,j,n,l,lam,bk,Np,Nq)
    !
    ! INPUTs:
    !        o) Np(opt), Nq(opt), w, j, n, l, lam: Quantum numbers
    !        o) bk: one character = b (bra) or k (ket)
    !
    ! OUTPUT:
    !        o) pretty_braket: character type
    !
    implicit none
    !
    integer, intent(in)         :: w,j,n,l,lam
    integer,optional            :: Np, Nq
    character(len=1),intent(in) :: bk
    character(len=250)          :: pretty_braket, aux
    !
    if ( (bk /= "b")  .and. (bk /= "k") ) &
         STOP "ERROR! pretty_braket: bk must be equal to b (bra) ot k (ket) !!! "
    !
    if ( bk == 'b' ) then
       write(pretty_braket,'(A2)') "< "
    else
       write(pretty_braket,'(A2)') "| "
    end if
    !
    if (present(Np)) then
       write(aux,'(I10)') Np
       pretty_braket=trim(pretty_braket) // " [Np=" // trim(adjustl(aux)) // "]"
    endif
    !
    write(aux,'(I10)') w
    pretty_braket=trim(pretty_braket) // " w=" // trim(adjustl(aux))
    !
    write(aux,'(I10)') j
    pretty_braket=trim(pretty_braket) // " j=" // trim(adjustl(aux)) // ";"
    !
    if (present(Nq)) then
       write(aux,'(I10)') Nq
       pretty_braket=trim(pretty_braket) // " [Np=" // trim(adjustl(aux)) // "]"
       aux=""
    endif
    !
    write(aux,'(I10)') n
    pretty_braket=trim(pretty_braket) // " n=" // trim(adjustl(aux))
    !
    write(aux,'(I10)') l
    pretty_braket=trim(pretty_braket) // " l=" // trim(adjustl(aux)) // ";"
    !
    write(aux,'(I10)') lam
    pretty_braket=trim(pretty_braket) // " lambda=" // trim(adjustl(aux))
    !
    if ( bk == 'b' ) then
       pretty_braket=trim(pretty_braket) // " |"
    else
       pretty_braket=trim(pretty_braket) // " >"
    endif
    !
  end function pretty_braket
  !
  !*****************************************************************************************
  !
  subroutine read_expdat(file_nm,data,dim,nt)
    !
    ! Read the input file.
    ! Line style:
    ! 0,0,2,0,0;  1,0,1,1,1;  3866.0;  0.47
    !
    implicit none
    !
    character(len=50), intent(in):: file_nm
    integer,intent(in):: dim,nt
    type(exp_point):: data(1:dim)
    integer:: line
    !
    open(unit=nt,file=trim(file_nm),status='old',action='read')
    !
    do line=1,dim
       !
       read(nt,123) data(line)%ist(1),data(line)%ist(2),data(line)%ist(3), &
            data(line)%ist(4),data(line)%ist(5), data(line)%fst(1),data(line)%fst(2), &
            data(line)%fst(3),data(line)%fst(4),data(line)%fst(5),data(line)%energy, &
            data(line)%intensity
       data(line)%ist(1)=Npval-2*data(line)%ist(1)
       data(line)%fst(1)=Npval-2*data(line)%fst(1)
       !
       if ( mod(data(line)%ist(2),2) == 0 ) then
          !
          data(line)%i_pos = find_pos( &
               (/ data(line)%ist(1),data(line)%ist(2),data(line)%ist(3), &
               data(line)%ist(4),data(line)%ist(5) /) , basis_para )
          !
          data(line)%f_pos = find_pos( &
               (/ data(line)%fst(1),data(line)%fst(2),data(line)%fst(3), &
               data(line)%fst(4),data(line)%fst(5) /) , basis_para )
          !
       else
          !
          data(line)%i_pos = find_pos( &
               (/ data(line)%ist(1),data(line)%ist(2),data(line)%ist(3), &
               data(line)%ist(4),data(line)%ist(5) /) , basis_ortho )
          !
          data(line)%f_pos = find_pos( &
               (/ data(line)%fst(1),data(line)%fst(2),data(line)%fst(3), &
               data(line)%fst(4),data(line)%fst(5) /) , basis_ortho )
          !
       endif
       !
    enddo
    !
123 format(10I2,F7.1,F4.2)
    !
    close(nt)
    !
  end subroutine read_expdat
  !
  !*****************************************************************************************
  !
  function exp_lines(file_nm,nt)
    !
    ! Return the number of lines of a file
    !
    implicit none
    !
    character(len=50),intent(in):: file_nm
    integer,intent(in):: nt
    integer:: exp_lines
    !
    open(unit=nt,file=trim(file_nm),status='old',action='read')
    !
    exp_lines=0
    do
       read(nt,*,end=12)
       exp_lines=exp_lines+1
    enddo
    !
12  close(nt)
    !
  end function exp_lines
  !
  !*****************************************************************************************
  !
  subroutine build_ham(H,bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1, basis)
    !
    !
    !
    implicit none
    !
    double precision, intent(in):: bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1
    integer,intent(in):: basis(:,:)
    integer:: ldim(1:2), udim(1:2)
    double precision,intent(inout):: H(:,:)
    integer:: i
    !
    ldim = lbound(basis)
    udim = ubound(basis)
    !
    if ( (ldim(1)/=1) .or. (udim(1)/=5) ) stop "ERROR! build_ham: 5 &
         &quantum numbers are needed !!! "
    !
    H=0.0d0
    !
    do i=ldim(2),udim(2)
       !
       H(i,i)=H(i,i)+bet*RME_Casimir_SOp4(basis(1,i),basis(2,i),basis(1,i),basis(2,i)) + &
            !
            gam* RME_Casimir_SOp3(basis(1,i),basis(2,i),basis(1,i),basis(2,i)) + &
            !
            gam2*(RME_Casimir_SOp3(basis(1,i),basis(2,i),basis(1,i),basis(2,i))**2.0d0) + &
            !
            kap*RME_Casimir_SOp4(basis(1,i),basis(2,i),basis(1,i),basis(2,i))* &
            RME_Casimir_SOp3(basis(1,i),basis(2,i),basis(1,i),basis(2,i)) + &
            !
            !
            a*RME_Casimir_uq3(basis(3,i),basis(4,i),basis(3,i),basis(4,i)) + &
            !
            b*(RME_Casimir_uq3(basis(3,i),basis(4,i),basis(3,i),basis(4,i))**2.0d0) + &
            !
            c*RME_Casimir_SOq3(basis(3,i),basis(4,i),basis(3,i),basis(4,i)) + &
            !
            !
            v1*RME_Casimir_uq3(basis(3,i),basis(4,i),basis(3,i),basis(4,i))* &
            RME_Casimir_SOp4(basis(1,i),basis(2,i),basis(1,i),basis(2,i))
       !
    enddo
    !
    if ( mod(basis(2,1),2) == 0 ) then! para case
       !
       H = H + d*SOq4(basis(5,1))%para + Qpq*QpQq(basis(5,1))%para + &
            Qpqw*QpQqW(basis(5,1))%para
       !
    else
       !
       H = H + d*SOq4(basis(5,1))%ortho + Qpq*QpQq(basis(5,1))%ortho + &
            QpQW*qpqqw(basis(5,1))%ortho
       !
    endif
    !
  end subroutine build_ham
  !
  !*****************************************************************************************
  !
  function find_pos(state,basis,iprint)
    !
    ! Find the position of /state/ in the /basis/ given
    !
    integer,intent(in):: state(1:5)
    integer,intent(in):: basis(:,:)
    integer,optional::iprint
    integer:: ldim(1:2), udim(1:2)
    integer:: find_pos
    !
    ldim=lbound(basis)
    udim=ubound(basis)
    !
    if ( (ldim(1)/=1) .or. (udim(1)/=5) ) stop "ERROR! build_ham: 5 &
         &quantum numbers are needed !!! "
    !
    find_pos = 1
    !
    do while ( (state(1)/=basis(1,find_pos)) .or. (state(2)/=basis(2,find_pos)) .or. &
         (state(3)/=basis(3,find_pos)) .or. (state(4)/=basis(4,find_pos)) .or. &
         (state(5)/=basis(5,find_pos)) .and. (find_pos.le.udim(2)) )
       !
       find_pos=find_pos+1
       !
    enddo
    if (present(iprint)) then
       !
       write(*,'(A,A,I8)')trim(pretty_braket(state(1),state(2),state(3), &
            state(4),state(5),'k')), " position: ", find_pos
       !
    endif
    !
  end function find_pos
  !
  !*****************************************************************************************
  !
  subroutine eigensystem(paraE,orthoE,lambda_max, &
       bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1,  &
       avec )
    !
    ! 26/02/2021
    ! avec = .False. (default)
    ! avec = .True. : H ---> Eigenfunctions
    !
    implicit none
    !
    !type(matrix),intent(inout):: Ham(1:lambda_max)
    double precision, intent(inout)::  paraE(1:sum(dim_para)), orthoE(1:sum(dim_ortho))
    integer,intent(in):: lambda_max
    double precision, intent(in):: bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1
    integer::i,j,k, state(1:5)
    double precision, allocatable:: aux(:)
    double precision:: ZPE
    !
    ! Optional argument to compute eigenvectors
    logical,optional:: avec
    ! Auxiliary memory to store temporally the coefficients
    !
    type(matrix):: auxV(0:lambda_max)
    !
#ifdef __INTEL_COMPILER
    double precision, allocatable:: vectors(:,:)
    call mkl_set_dynamic(0)
#endif
    ! Default value of avec...
    if (.not. present(avec)) then
       avec = .false.
    elseif (avec .eqv. .true.) then
       ! Allocate auxV:
       do i=0,lambda_max
          !
          allocate( AuxV(i)%para(1:dim_para(i),1:dim_para(i)), &
               AuxV(i)%ortho(1:dim_ortho(i),1:dim_ortho(i)) )
          !
          AuxV(i)%para = 0.0d0
          AuxV(i)%ortho= 0.0d0
          !
       enddo
       !
    endif
    !
#ifdef __GFORTRAN__
    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(aux,j,k,state)
#elif defined __INTEL_COMPILER
    !$OMP PARALLEL DO DEFAULT(SHARED), PRIVATE(aux,j,k,state,vectors)
#endif
    do i=0,lambda_max
       !
       ! Build the Hamiltonians
       call build_ham(Ham(i)%para,bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1,basis_para( &
            1:5,ijk_para(i):ijk_para(i)+dim_para(i)-1) )
       !
       call build_ham(Ham(i)%ortho,bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1,basis_ortho( &
            1:5,ijk_ortho(i):ijk_ortho(i)+dim_ortho(i)-1) )
       !
       ! Diagonalize - para
       allocate(aux(1:dim_para(i)))
       !
#ifdef __GFORTRAN__
       call la_syevr(A=Ham(i)%para,W=aux,JOBZ='V',UPLO='U')
#elif defined __INTEL_COMPILER
       allocate(vectors(1:dim_para(i),1:dim_para(i)))
       call syevr(A=Ham(i)%para,W=aux,UPLO='U',Z=vectors)
       Ham(i)%para=vectors
       deallocate(vectors)
#endif
       do j =1,dim_para(i)
          ! (:,fix_col)
          state = assig_state(Ham(i)%para(:,j),basis_para( &
               1:5,ijk_para(i):ijk_para(i)+dim_para(i)-1) )
          k = find_pos(state,basis_para)
          paraE(k) = aux(j)
          ! Ordering eigenvectors if avec=.true.
          if ( avec .eqv. .true.) then
             auxV(i)%para(:,k-sum(dim_para(0:i-1))) = Ham(i)%para(:,j)
          endif
       enddo
       deallocate(aux)
       !
       ! Diagonalize - ortho
       allocate(aux(1:dim_ortho(i)))
       !
#ifdef __GFORTRAN__
       call la_syevr(A=Ham(i)%ortho,W=aux,JOBZ='V',UPLO='U')
#elif defined __INTEL_COMPILER
       allocate(vectors(1:dim_ortho(i),1:dim_ortho(i)))
       call syevr(A=Ham(i)%ortho,W=aux,UPLO='U',Z=vectors)
       Ham(i)%ortho=vectors
       deallocate(vectors)
#endif
       do j =1,dim_ortho(i)
          ! (:,fix_col)
          state = assig_state(Ham(i)%ortho(:,j),basis_ortho( &
               1:5,ijk_ortho(i):ijk_ortho(i)+dim_ortho(i)-1) )
          k = find_pos(state,basis_ortho)
          orthoE(k) = aux(j)
          ! Ordering eigenvectors if avec=.true.
          if ( avec .eqv. .true.) then
             auxV(i)%ortho(:,k-sum(dim_ortho(0:i-1))) = Ham(i)%ortho(:,j)
          endif
       enddo
       deallocate(aux)
       !
       ! If avec ---> eigenvectors are stored in Hmatrix (so ... the
       ! Hamiltonian is killed ...)
       if (avec .eqv. .true.) then
          Ham(i)%para = auxV(i)%para
          Ham(i)%ortho= auxV(i)%ortho
       endif
       !
    enddo
    !
    ZPE = minval(paraE)
    paraE = paraE - ZPE
    orthoE = orthoE - ZPE
    !
  end subroutine eigensystem    
  !
  !*****************************************************************************************
  !
  function assig_state(WF,basis)
    !
    ! Give the WF components respect a basis and a state is asignated
    !
    implicit none
    !
    double precision, intent(in)::WF(:)
    integer,intent(in)::basis(:,:)
    integer:: assig_state(1:5), pos(1)
    !
    pos = maxloc( WF*WF ) ! position of the maximun square component
    assig_state(:) = basis(:,pos(1))
    !
  end function assig_state
  !
  !*****************************************************************************************
  !
  function chi2(params)
    !
    ! params(1) ---> bet
    ! params(2) ---> gam
    ! params(3) ---> gam2
    ! params(4) ---> kap
    !
    ! params(5) ---> a
    ! params(6) ---> b
    ! params(7) ---> c
    ! params(8) ---> d
    !
    ! params(9) ---> Qpq
    ! params(10) --> QpqW
    ! params(11) --> v1
    !
    implicit none
    !
    double precision,intent(in):: params(1:11)
    double precision:: chi2,aux
    double precision:: paraE(1:sum(dim_para)), orthoE(1:sum(dim_ortho))
    integer:: i
    integer::miflag
    common/minuit_iflag/miflag
    !
    chi2 = 0.0d0
    !
    call eigensystem( paraE,orthoE,lambda_max, &
         params(1), params(2), params(3), params(4), &
         params(5), params(6), params(7), params(8), &
         params(9), params(10), params(11) )
    !
    ! Counpute chi2
    do i=1,total_exp
       !
       aux=0.0d0
       !
       if ( mod(exp_data(i)%ist(2),2) == 0 ) then ! para-state
          !
          aux = exp_data(i)%energy - (paraE(exp_data(i)%f_pos) - paraE(exp_data(i)%i_pos))
          !
       else !ortho-state
          !
          aux = exp_data(i)%energy - (orthoE(exp_data(i)%f_pos) - orthoE(exp_data(i)%i_pos))
          !
       endif
       !
       chi2 = chi2 + aux*aux
       !
    enddo
    !
    if (miflag == 3) then
       !
       write(*,'(/,A/)') "Residuals: exp - calc"
       do i=1,total_exp
          !
          if ( mod(exp_data(i)%ist(2),2) == 0 ) then ! para-state
             !
             write(*,'(A,F7.3)') " Para:"// trim(pretty_braket(exp_data(i)%ist(1), &
                  exp_data(i)%ist(2),exp_data(i)%ist(3),exp_data(i)%ist(4), &
                  exp_data(i)%ist(5),'k'))//' ---> '// trim(pretty_braket( &
                  exp_data(i)%fst(1),exp_data(i)%fst(2),exp_data(i)%fst(3), &
                  exp_data(i)%fst(4),exp_data(i)%fst(5),'k'))// ' : ', &
                  exp_data(i)%energy - (paraE(exp_data(i)%f_pos) - &
                  paraE(exp_data(i)%i_pos))
             !
          else !ortho-state
             !
             write(*,'(A,F7.3)') "Ortho:"// trim(pretty_braket(exp_data(i)%ist(1), &
                  exp_data(i)%ist(2),exp_data(i)%ist(3),exp_data(i)%ist(4), &
                  exp_data(i)%ist(5),'k')) //' ---> '// trim(pretty_braket( &
                  exp_data(i)%fst(1),exp_data(i)%fst(2),exp_data(i)%fst(3), &
                  exp_data(i)%fst(4),exp_data(i)%fst(5),'k'))// ' : ', &
                  exp_data(i)%energy - (orthoE(exp_data(i)%f_pos) - &
                  orthoE(exp_data(i)%i_pos))
             !
          endif
          !
       enddo
       !
       write(*,*)
       write(*,*) "    Beta:", params(1)
       write(*,*) "   Gamma:", params(2)
       write(*,*) "  Gamma2:", params(3)
       write(*,*) "   Kappa:", params(4)
       write(*,*)
       write(*,*) "       a:", params(5)
       write(*,*) "       b:", params(6)
       write(*,*) "       c:", params(7)
       write(*,*) "       d:", params(8)
       write(*,*)
       write(*,*) "     Qpq:", params(9)
       write(*,*) "    QpqW:", params(10)
       write(*,*) "      v1:", params(11)
       
       !
    endif
    !
  end function chi2
  !
  !*****************************************************************************************
  !
  subroutine FCN(npar,grad,fval,xval,iflag,chi2)
    !
    ! Minuit function
    !
    implicit none
    !
    double precision:: grad(*), xval(*), fval
    integer:: iflag, npar
    double precision, external:: chi2
    integer:: miflag
    common/minuit_iflag/miflag
    miflag=iflag
    !
    fval=chi2(xval)
    !
    write(*,31) fval, sqrt(fval/dble(total_exp-npar)), sqrt(fval/dble(total_exp))
31  format('CHI2 =',F17.2,', RMS = ', F14.4, ', SIGMA = ', F14.4)
    !
  end subroutine FCN
  !
  !*****************************************************************************************
  !
  function RME_np_x_Dq_1(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    ! < state 1 || [ np x Dq' ]^(1) || state 2 >
    !
    integer, intent(in):: w1,J1,n1,L1,lam1, w2,J2,n2,L2,lam2
    double precision:: RME_np_x_Dq_1
    !
    RME_np_x_Dq_1 = sqrt( 3.0d0 * dble(2*lam1+1)*dble(2*lam2+1) ) * &
         wigner_9j( j1,l1,lam1, j2,l2,lam2, 0,1,1 ) * &
         RME_Dq_prima(n1,l1,n2,l2) * RME_np(w1,j1,w2,j2)
    !
  end function RME_np_x_Dq_1
  !
  !*****************************************************************************************
  !
  function RME_Qp2_x_Dq_1(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    !  < state 1 || [ Qp2 x Dq' ]^(1) || state 2 >
    !
    integer, intent(in):: w1,J1,n1,L1,lam1, w2,J2,n2,L2,lam2
    double precision:: RME_Qp2_x_Dq_1
    !
    RME_Qp2_x_Dq_1 = sqrt( 3.0d0 * dble(2*lam1+1)*dble(2*lam2+1) ) * &
         wigner_9j( j1,l1,lam1, j2,l2,lam2, 2,1,1 ) * &
         RME_Dq_prima(n1,l1,n2,l2) * RME_Qp2(w1,j1,w2,j2)
    !
  end function RME_Qp2_x_Dq_1
  !
  !*****************************************************************************************
  !
  function RME_Qp2_x_nq_2(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    !  < state 1 || [ Qp2 x nq ]^(2) || state 2 > ! Quadruple 
    !
    integer, intent(in):: w1,J1,n1,L1,lam1, w2,J2,n2,L2,lam2
    double precision:: RME_Qp2_x_nq_2
    !
    RME_Qp2_x_nq_2 = sqrt( 5.0d0 * dble(2*lam1+1)*dble(2*lam2+1) ) * &
         wigner_9j(j1,l1,lam1, j2,l2,lam2, 2,0,2) * &
         RME_Casimir_Uq3(n1,l1,n2,l2) * RME_Qp2(w1,j1,w2,j2)
    !
  end function RME_Qp2_x_nq_2
  !
  !*****************************************************************************************
  !
  subroutine StoreEigenMatel(intop)
    !
    ! Computes < psi1,lambda1, S | T^(t)| psi2, lambda2, S >
    ! and store them in intop(len(exp_Data),1:3)
    !
    double precision, intent(out):: intop(:,:)
    integer:: &
         trans, & ! loop over all experimental transition
         i_pos, f_pos
    !
    do trans=1,total_exp
       !
       i_pos=exp_data(trans)%i_pos
       f_pos=exp_data(trans)%f_pos
       ! EigenExpected(coefs1,lamb1,coefs2,lamb2,S,RME_fun)
       if ( mod(exp_data(trans)%ist(2),2) == 0 ) then! para-case
          !
          intop(trans,1) = EigenExpected( &
               Ham( exp_data(trans)%ist(5) )%para(:,i_pos- &
               sum(dim_para(0:exp_data(trans)%ist(5)-1))), & ! coefs1
               exp_data(trans)%ist(5),                     & ! lamb1
               Ham( exp_data(trans)%fst(5) )%para(:,f_pos- &
               sum(dim_para(0:exp_data(trans)%fst(5)-1))), & ! coef2
               exp_data(trans)%fst(5),                     & ! lamb2
               0,  RME_np_x_Dq_1)
          !
          intop(trans,2) = EigenExpected( &
               Ham( exp_data(trans)%ist(5) )%para(:,i_pos- &
               sum(dim_para(0:exp_data(trans)%ist(5)-1))), & ! coefs1
               exp_data(trans)%ist(5),                     & ! lamb1
               Ham( exp_data(trans)%fst(5) )%para(:,f_pos- &
               sum(dim_para(0:exp_data(trans)%fst(5)-1))), & ! coef2
               exp_data(trans)%fst(5),                     & ! lamb2
               0,  RME_Qp2_x_Dq_1)
          !
          intop(trans,3) = EigenExpected( &
               Ham( exp_data(trans)%ist(5) )%para(:,i_pos- &
               sum(dim_para(0:exp_data(trans)%ist(5)-1))), & ! coefs1
               exp_data(trans)%ist(5),                     & ! lamb1
               Ham( exp_data(trans)%fst(5) )%para(:,f_pos- &
               sum(dim_para(0:exp_data(trans)%fst(5)-1))), & ! coef2
               exp_data(trans)%fst(5),                     & ! lamb2
               0,  RME_Qp2_x_nq_2)
          !
       else ! ortho-case
          !
          intop(trans,1) = EigenExpected( &
               Ham( exp_data(trans)%ist(5) )%ortho(:,i_pos- &
               sum(dim_ortho(0:exp_data(trans)%ist(5)-1))), & ! coefs1
               exp_data(trans)%ist(5),                      & ! lamb1
               Ham( exp_data(trans)%fst(5) )%ortho(:,f_pos- &
               sum(dim_ortho(0:exp_data(trans)%fst(5)-1))), & ! coef2
               exp_data(trans)%fst(5),                      & ! lamb2
               1,  RME_np_x_Dq_1) 
          !
          intop(trans,2) = EigenExpected( &
               Ham( exp_data(trans)%ist(5) )%ortho(:,i_pos- &
               sum(dim_ortho(0:exp_data(trans)%ist(5)-1))), & ! coefs1
               exp_data(trans)%ist(5),                      & ! lamb1
               Ham( exp_data(trans)%fst(5) )%ortho(:,f_pos- &
               sum(dim_ortho(0:exp_data(trans)%fst(5)-1))), & ! coef2
               exp_data(trans)%fst(5),                      & ! lamb2
               1,  RME_Qp2_x_Dq_1)
          !
          intop(trans,3) = EigenExpected( &
               Ham( exp_data(trans)%ist(5) )%ortho(:,i_pos- &
               sum(dim_ortho(0:exp_data(trans)%ist(5)-1))), & ! coefs1
               exp_data(trans)%ist(5),                      & ! lamb1
               Ham( exp_data(trans)%fst(5) )%ortho(:,f_pos- &
               sum(dim_ortho(0:exp_data(trans)%fst(5)-1))), & ! coef2
               exp_data(trans)%fst(5),                      & ! lamb2
               1,  RME_Qp2_x_nq_2)
          !
       endif
       !
    enddo
    !
  end subroutine StoreEigenMatel
  !
  !*****************************************************************************************
  !
  function EigenExpected(coefs1,lamb1,coefs2,lamb2,S,RME_fun)
    !
    ! < psi1,lambda1, S | T^(t)| psi2, lambda2, S >
    ! Input:
    ! o) Coefs1: length = len(Block_lamb1)
    ! o) lamb1
    ! o) Coefs2: length = len(Block_lamb2)
    ! o) lamb2
    ! o) S : 0 if para, 1 if ortho
    ! o) RME_fun (w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    ! Output:
    ! o) < psi1,lambda1, S | T^(t)| psi2, lambda2, S >
    !
    !
    double precision, intent(in) :: coefs1(:), coefs2(:)
    integer, intent(in)         :: lamb1, lamb2, S
    double precision, external  :: RME_fun
    double precision            :: EigenExpected
    integer                     :: i,j
    !
    EigenExpected = 0.0d0
    !
    if (S==0) then
       !
       do i = 1,dim_para(lamb1) !ijk_para(lamb1),ijk_para(lamb1)+dim_para(lamb1)
          do j = 1,dim_para(lamb2) !ijk_para(lamb2),ijk_para(lamb2)+dim_para(lamb2)
             !
             EigenExpected = EigenExpected +   &
                  coefs1(i)* &
                  coefs2(j)* &
                  RME_fun( &
                  basis_para(1,i+ijk_para(lamb1)-1), &
                  basis_para(2,i+ijk_para(lamb1)-1), &
                  basis_para(3,i+ijk_para(lamb1)-1), &
                  basis_para(4,i+ijk_para(lamb1)-1), &
                  basis_para(5,i+ijk_para(lamb1)-1), &
                  basis_para(1,j+ijk_para(lamb2)-1), &
                  basis_para(2,j+ijk_para(lamb2)-1), &
                  basis_para(3,j+ijk_para(lamb2)-1), &
                  basis_para(4,j+ijk_para(lamb2)-1), &
                  basis_para(5,j+ijk_para(lamb2)-1))
             !
          enddo
       enddo
       !
    else if (S==1) then
       !
       do i = 1,dim_ortho(lamb1) !ijk_ortho(lamb1),ijk_ortho(lamb1)+dim_ortho(lamb1)
          do j = 1,dim_ortho(lamb2) !ijk_ortho(lamb2),ijk_ortho(lamb2)+dim_ortho(lamb2)
             !
             EigenExpected = EigenExpected + &
                  coefs1(i)* &
                  coefs2(j)* &
                  RME_fun( &
                  basis_ortho(1,i+ijk_ortho(lamb1)-1), &
                  basis_ortho(2,i+ijk_ortho(lamb1)-1), &
                  basis_ortho(3,i+ijk_ortho(lamb1)-1), &
                  basis_ortho(4,i+ijk_ortho(lamb1)-1), &
                  basis_ortho(5,i+ijk_ortho(lamb1)-1), &
                  basis_ortho(1,j+ijk_ortho(lamb2)-1), &
                  basis_ortho(2,j+ijk_ortho(lamb2)-1), &
                  basis_ortho(3,j+ijk_ortho(lamb2)-1), &
                  basis_ortho(4,j+ijk_ortho(lamb2)-1), &
                  basis_ortho(5,j+ijk_ortho(lamb2)-1))
             !
          enddo
       enddo
       !
    end if
    !
  end function EigenExpected
  !
  !*****************************************************************************************
  !
  function chi2_int(params)
    !
    ! params(1) ---> npo
    ! params(2) ---> dnD
    ! params(3) ---> dQD
    ! params(4) ---> qQn
    !
    implicit none
    !
    double precision,intent(in):: params(1:4)
    double precision:: chi2_int, comp_int
    integer:: i
    integer::miflag
    common/minuit_iflag/miflag
    !
    chi2_int = 0.0d0
    !
    do i=1,total_exp
       !

       !
       chi2_int = chi2_int + ( exp_data(i)%intensity - &
            compute_transition(params,exp_data(i),i) )**2.0d0
       !
    enddo
    !
    if (miflag == 3) then
       !
       write(*,'(/,A/)') "Residuals: exp - calc"
       do i=1,total_exp
          !
          comp_int = compute_transition(params,exp_data(i),i)
          !
          if ( mod(exp_data(i)%ist(2),2) == 0 ) then ! para-state
             !
             write(*,'(A,F8.4,A,F8.4)') " Para:"// trim(pretty_braket(exp_data(i)%ist(1), &
                  exp_data(i)%ist(2),exp_data(i)%ist(3),exp_data(i)%ist(4), &
                  exp_data(i)%ist(5),'k'))//' ---> '// trim(pretty_braket( &
                  exp_data(i)%fst(1),exp_data(i)%fst(2),exp_data(i)%fst(3), &
                  exp_data(i)%fst(4),exp_data(i)%fst(5),'k'))// ' : ', &
                  (exp_data(i)%intensity - comp_int)*100.0/exp_data(i)%intensity, &
                  "% ---", comp_int 
             !
          else !ortho-state
             !
             write(*,'(A,F8.4,A,F8.4)') "Ortho:"// trim(pretty_braket(exp_data(i)%ist(1), &
                  exp_data(i)%ist(2),exp_data(i)%ist(3),exp_data(i)%ist(4), &
                  exp_data(i)%ist(5),'k')) //' ---> '// trim(pretty_braket( &
                  exp_data(i)%fst(1),exp_data(i)%fst(2),exp_data(i)%fst(3), &
                  exp_data(i)%fst(4),exp_data(i)%fst(5),'k'))// ' : ', &
                  (exp_data(i)%intensity - comp_int)*100.d0/exp_data(i)%intensity, &
                  "% ---", comp_int 
          endif
          !
       enddo
       !
       write(*,*)
       write(*,*) "    npo npara/northo  :", params(1)
       write(*,*) "    dnD [np x D'q]^(1):", params(2)
       write(*,*) "    dQD [Qp x D'2]^(1):", params(3)
       write(*,*) "    qQn [Qp x nq ]^(2):", params(4)       
       !
    endif
    !
  end function chi2_int
  !  
  !*****************************************************************************************
  !
  function compute_ProbTrans(npo,transition)
    !
    ! pi = nk (2lamb_i+1) exp(-Ei/KbT) / Sum_j[ 2lamb_j+1) exp(-Ej/KbT)]
    !
    ! Inputs:
    ! o) npo: ratio between para and ortho
    ! o) transitions: exp_point type
    !
    double precision, intent(in):: npo
    type(exp_point),intent(in):: transition
    double precision:: compute_ProbTrans, suma
    integer :: lamb, i
    !
    suma = 0.0d0
    if (mod(transition%ist(2),2)==0) then !para-case
       !
       do lamb=0,lambda_max
          do i = ijk_para(lamb),ijk_para(lamb)+dim_para(lamb)-1
             !
             suma = suma + dble(2*lamb+1)*exp(-EnergiesPara(i) / temp)
             !
          enddo
       enddo
       !
       compute_ProbTrans = npo * dble(2*transition%ist(5)+1) * &
            exp(-EnergiesPara(transition%i_pos) / Temp )
       !
    else !ortho-case
       !
       do lamb=0,lambda_max
          do i = ijk_ortho(lamb),ijk_ortho(lamb)+dim_ortho(lamb)-1
             !
             suma = suma + dble(2*lamb+1)*exp(-EnergiesOrtho(i) / temp)
             !
          enddo
       enddo
       !
       compute_ProbTrans = (1.0d0/npo) * dble(2*transition%ist(5)+1) * &
            exp(-EnergiesOrtho(transition%i_pos) / Temp )
       !
    endif
    !
  end function compute_ProbTrans
  !  
  !*****************************************************************************************
  !
  function compute_transition(params,transition,i)
    !
    ! Full transition:
    ! Inputs:
    ! o) params:
    !    params(1) ---> npo
    !    params(2) ---> dnD
    !    params(3) ---> dQD
    !    params(4) ---> qQn
    ! o) transition: type exp_data
    ! o) i: position in exp_data array
    !
    !
    double precision,intent(in):: params(1:4)
    type(exp_point), intent(in):: transition
    integer, intent(in)        :: i
    double precision           :: compute_transition
    !
    if( mod(transition%ist(2),2)==0 ) then !para-case
       !
       compute_transition = compute_ProbTrans(params(1),transition)* &
            ( EnergiesPara( transition%f_pos ) - EnergiesPara( transition%i_pos ) ) * &
            abs(params(2)*intop(i,1)+params(3)*intop(i,2)+params(4)*intop(i,3))
       !
    else !ortho-case
       !
       compute_transition = compute_ProbTrans(params(1),transition)* &
            ( EnergiesOrtho( transition%f_pos ) - EnergiesOrtho( transition%i_pos ) ) * &
            abs(params(2)*intop(i,1)+params(3)*intop(i,2)+params(4)*intop(i,3))
       !
    endif
  end function compute_transition
  !  
  !*****************************************************************************************
  !
  subroutine FCN_int(npar,grad,fval,xval,iflag,chi2)
    !
    ! Minuit function
    !
    implicit none
    !
    double precision:: grad(*), xval(*), fval
    integer:: iflag, npar
    double precision, external:: chi2
    integer:: miflag
    common/minuit_iflag/miflag
    miflag=iflag
    !
    fval=chi2(xval)
    !
    write(*,31) fval, sqrt(fval/dble(total_exp-npar)), sqrt(fval/dble(total_exp))
31  format('CHI2 =',F17.2,', RMS = ', F14.4, ', SIGMA = ', F14.4)
    !
  end subroutine FCN_INT
  
  !  
  !*****************************************************************************************
  !
end module MOD_Up_x_Uq




