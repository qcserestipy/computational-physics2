module cp_mod
    implicit none
    save

!##########################
!#    type definitions    #
!##########################

    type :: point
        real :: x = 0.0d0
        real :: y = 0.0d0
        real :: z = 0.0d0
    end type point

    type :: atom
        character(len=2) :: name ='H'
        integer :: number
        real :: mass = 1.0d0
        real :: charge = 0.0d0
        type( point ) :: coords
    end type atom

    type :: system
        type( atom ) :: atoms
    end type system


!##############################
!#    function definitions    #
!##############################
    contains
        subroutine readxyzfile(filename, listofatoms)
            implicit none

            type( atom ), dimension(:), allocatable, intent(inout) :: listofatoms
            character(len=72), intent(in) :: filename
            integer :: i, io_error, natoms

            open(unit=1,file=filename, iostat=io_error)
            if (io_error == 0) then
                read(1,*) natoms
                allocate(listofatoms(natoms))
                read(1,*) 

                do i = 1, natoms
                    read(1,*) listofatoms(i)%name, &
                    listofatoms(i)%coords%x, &
                    listofatoms(i)%coords%y, &
                    listofatoms(i)%coords%z
                enddo

                write(*,*) "Reading file, getting atoms with following coordinates:" 
                write(*,*) "-------------------------------------------------------------"
                do i = 1, natoms
                    write(*,*) "|", listofatoms(i)%name, &
                    "|", listofatoms(i)%coords%x, &
                    "|", listofatoms(i)%coords%y, &
                    "|", listofatoms(i)%coords%z, "|"
                enddo
                write(*,*) "-------------------------------------------------------------"
            else
                write(*,*) 'Opening the file caused an error', io_error
            end if
            close(1)
        end subroutine readxyzfile

        subroutine writexyzfile(filename, listofatoms,length)
            implicit none

            type( atom ), dimension(:), allocatable, intent(inout) :: listofatoms
            character(len=72), intent(in) :: filename
            integer, intent(in) :: length
            integer :: i, io_error
            
            open(unit=1,file=filename, iostat=io_error)
                write(1,*) length
                write(1,*)
                do i = 1, length
                    write(1,*) listofatoms(i)%name, &
                    listofatoms(i)%coords%x, &
                    listofatoms(i)%coords%y, &
                    listofatoms(i)%coords%z
                enddo
            close(1)
        end subroutine writexyzfile

        subroutine printmatrix(a)
            implicit none
            integer :: i,j
            real, dimension(:,:) :: a
            write(*,*)            
            do i = lbound(a,1), ubound(a,1)
               write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
            end do
         end subroutine printmatrix

end module cp_mod