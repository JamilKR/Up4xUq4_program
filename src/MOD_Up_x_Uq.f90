module MOD_Up_x_Uq
  !
  use MOD_matfun
  use MOD_Up4
  use MOD_Uq4
  !
  implicit none
  !
  integer:: total_dim
  integer:: para_dim
  integer:: ortho_dim
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
    !         o) basis_para (1:5,0:lambda_max): basis  para-blocks
    !         o) basis_ortho(0:5,0:lambda_max): basis ortho-blocks
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
         &5 quantum number are needed !!! "
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
end module MOD_Up_x_Uq




