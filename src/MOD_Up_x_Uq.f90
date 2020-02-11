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
  function build_Qp_x_Qq_0_matrix(basis,dim,iprint)
    !
    ! sqrt(5) factor omited!
    !
    implicit none
    !
    integer,intent(in):: dim
    integer,intent(in):: basis(1:5,1:dim)
    integer,optional:: iprint
    double precision::   build_Qp_x_Qq_0_matrix(1:dim,1:dim)
    integer:: i,j
    !
    build_Qp_x_Qq_0_matrix=0.0d0
    !
    do i=1,dim ! bra-index
       !
       if(present(iprint)) write(*,"(7(A,I3),A)") "<[",Npval,"],",basis(1,i),",",&
            basis(2,i),";[",Nqval,"],",basis(3,i),",",basis(4,i),";",basis(5,i),"|"
       !
       do j=1,dim ! ket-index
          !
          if(present(iprint)) write(*,"(T34,7(A,I3),A)") "|[",Npval,"],",&
               basis(1,j),",",basis(2,j),";[",Nqval,"],",basis(3,j),",",basis(4,j),";", &
               basis(5,i),">"
          !
          if ( basis(5,i) /= basis(5,j) ) cycle
          !
          build_Qp_x_Qq_0_matrix(i,j) = dble( (-1)**(basis(4,i)+basis(5,i)+basis(2,j)) ) * &
               wigner_6j(basis(2,i),basis(4,i),basis(5,i), basis(4,j),basis(2,j),2) * &
               RME_Qq2(basis(3,i),basis(4,i),basis(3,j),basis(4,j)) * &
               RME_Qp2(basis(1,i),basis(2,i),basis(1,j),basis(2,j))
          !
          if(present(iprint))  write(*,"(T64,F15.4)") build_Qp_x_Qq_0_matrix(i,j)
          !
       end do
       !
    enddo      
    !
  end function build_Qp_x_Qq_0_matrix
  !
  !*****************************************************************************************
  !
end module MOD_Up_x_Uq




