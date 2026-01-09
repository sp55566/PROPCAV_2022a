! ----------------------------------------------------------------------------- !
!                       CREATED BY YIRAN SU 10/04/2017                          !
!             Linked List to Store PROPCAV-MESH Correlation!                    !
! ----------------------------------------------------------------------------- !

      module M_LNK_LST_P2M
        implicit none

        public  :: insertCor, getCorr  !, delCorr

        type :: COR_ENTRY
          integer :: ixp, irp
          integer :: ixm, irm
          double precision :: weight
          type(COR_ENTRY),pointer :: next
        end type COR_ENTRY

        type(COR_ENTRY),pointer,private :: first => null()
        type(COR_ENTRY),pointer,private :: last => null()

      contains

        ! Insert new items to the linked list
        subroutine insertCor(ixp,irp,ixm,irm,weight)
          implicit none
          type(COR_ENTRY),pointer :: new
          integer,intent(in) :: ixp, irp, ixm, irm
          double precision,intent(in) :: weight
          allocate(new)
          new%ixp = ixp
          new%irp = irp
          new%ixm = ixm
          new%irm = irm
          new%weight = weight
          new%next => null()
          if (.not. associated(first)) then
            first => new
            last => new
          else
            last%next => new
            last => new
          end if
        end subroutine

        ! Get the first item (pointer)
        function getCorr()
          implicit none
          type(COR_ENTRY),pointer :: getCorr
          getCorr => first
        end function

*        subroutine delCorr()
*          implicit none
*          type(COR_ENTRY),pointer :: temp
*          do while(associated(first))
*            temp => first%next
*            deallocate(first)
*            first => temp
*          end do
*          first => null()
*          last => null()
*        end subroutine

      end module


*      program test
*        use M_COR_P2M
*        implicit none
*        type(COR_ENTRY),pointer :: cor_p2m, temp
*        integer :: i
*        do i=1,11
*          call insertCor(i,i,i,i,10.0*i)
*        end do
*        cor_p2m => getCorr()
*        do while(associated(cor_p2m))
*          write(*,*) cor_p2m%ixp, cor_p2m%weight
*          temp => cor_p2m%next
*          cor_p2m => temp
*        end do
*      end program

