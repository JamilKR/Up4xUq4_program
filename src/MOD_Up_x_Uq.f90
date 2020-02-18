module MOD_Up_x_Uq
  !
  use MOD_matfun
  use MOD_Up4
  use MOD_Uq4
  USE F95_LAPACK, ONLY: LA_SYEVR
  !
  implicit none
  !
  integer, allocatable:: basis_para(:,:),basis_ortho(:,:), dim_para(:), dim_ortho(:), &
       ijk_para(:), ijk_ortho(:)
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
       QpQqW(:)    ! Qp x Qq SOp(4) matrices para and ortho
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
       read(nt,*) data(line)%ist(1),data(line)%ist(2),data(line)%ist(3), &
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
  function RME_sop4(w1,j1,n1,l1,lam1,w2,j2,n2,l2,lam2)
    !
    ! RME_sop4 = diagonal !!!
    ! RME_sop4 = w(w+2)
    !
    implicit none
    !
    integer,intent(in):: w1,j1,n1,l1,lam1, w2,j2,n2,l2,lam2
    double precision  :: RME_sop4
    !
    if ( (w1==w2) .and. (j1==j2) .and. (n1==n2) .and. (l1==l2) .and. (lam1==lam2) ) then
       !
       RME_sop4 = dble( w1*(w1+2) )
       !
    else
       !
       RME_sop4 = 0.0d0
       !
    endif
    !
  end function RME_sop4
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
    h=0.0d0
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
       h = h + d*SOq4(basis(5,1))%para + Qpq*QpQq(basis(5,1))%para + &
            Qpqw*QpQqW(basis(5,1))%para
       !
    else
       !
       h = h + d*SOq4(basis(5,1))%ortho + Qpq*QpQq(basis(5,1))%ortho + &
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
    integer:: ldim(1:2), udim(1:2), i
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
  subroutine eigensystem(Ham,paraE,orthoE,lambda_max, &
       bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1)
    !
    !
    !
    implicit none
    !
    type(matrix),intent(inout):: Ham(1:lambda_max)
    double precision, intent(inout)::  paraE(*), orthoE(*)
    integer,intent(in):: lambda_max
    double precision, intent(in):: bet,gam,gam2,kap, a,b,c,d, Qpq, Qpqw, v1
    integer::i,j,k, state(1:5)
    double precision, allocatable:: aux(:)
    !
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
       call la_syevr(A=Ham(i)%para,W=aux,JOBZ='V',UPLO='U')
       do j =ijk_para(i),ijk_para(i)+dim_para(i)-1
          ! (:,fix_col)
          state = assig_state(Ham(i)%para(:,j),basis_para( &
               1:5,ijk_para(i):ijk_para(i)+dim_para(i)-1) )
          k = find_pos(state,basis_para)
          paraE(k) = aux(j)
       enddo
       deallocate(aux)
       !
       ! Diagonalize - ortho
       allocate(aux(1:dim_ortho(i)))
       call la_syevr(A=Ham(i)%ortho,W=aux,JOBZ='V',UPLO='U')
       do j =ijk_ortho(i),ijk_ortho(i)+dim_ortho(i)-1
          ! (:,fix_col)
          state = assig_state(Ham(i)%ortho(:,j),basis_ortho( &
               1:5,ijk_ortho(i):ijk_ortho(i)+dim_ortho(i)-1) )
          k = find_pos(state,basis_ortho)
          orthoE(k) = aux(j)
       enddo
       deallocate(aux)
       !
       !
    enddo
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
    !
    !
    implicit none
    !
    double precision:: chi2
    !
    chi2 = 0.0d0
    !
  end function chi2
  !
  !*****************************************************************************************
  !
end module MOD_Up_x_Uq




