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

    type :: velocity
        real :: v_x = 0.0d0
        real :: v_y = 0.0d0
        real :: v_z = 0.0d0
    end type velocity

    type :: atom
        character(len=2) :: name ='H'
        integer :: number = 1
        real :: mass = 1.0d0
        real :: charge = 0.0d0
        type( point ) :: coords
        type( velocity ) :: velo
        real :: s_velo
    end type atom

    type :: system
        type( atom ), dimension(:), allocatable :: subsystem
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




        subroutine readtrjfile(filename, listofFrames,nframes)
            implicit none
            type( system ), dimension(:), allocatable, intent(inout) :: listofFrames
            type( atom ), dimension(:), allocatable :: listofatoms
            character(len=72), intent(in) :: filename
            integer, intent(in) :: nframes
            integer :: natoms, i, j, io_error

            allocate(listofFrames(nframes))

            
            open(unit=1,file=filename, iostat=io_error)
            if (io_error == 0) then
                do j = 1, nframes
                    read(1,*) natoms
                    allocate(listofatoms(natoms))
                    read(1,*) 

                    do i = 1, natoms
                        read(1,*) listofatoms(i)%name, &
                        listofatoms(i)%coords%x, &
                        listofatoms(i)%coords%y, &
                        listofatoms(i)%coords%z
                    enddo

                    write(*,*) "Frame: ", j 
                    write(*,*) "-------------------------------------------------------------"
                    do i = 1, natoms
                        write(*,*) "|", listofatoms(i)%name, &
                        "|", listofatoms(i)%coords%x, &
                        "|", listofatoms(i)%coords%y, &
                        "|", listofatoms(i)%coords%z, "|"
                    enddo
                    write(*,*) "-------------------------------------------------------------"
                    
                    listofFrames(j)%subsystem = listofatoms
                    deallocate(listofatoms)
                enddo
            else
                write(*,*) 'Opening the file caused an error', io_error
            end if
            close(1)
        end subroutine readtrjfile




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




        subroutine writetrjfile(filename, listofatoms,length)
            implicit none

            type( atom ), dimension(:), allocatable, intent(inout) :: listofatoms
            character(len=72), intent(in) :: filename
            integer, intent(in) :: length
            integer :: i, io_error
            
            open(unit=1,file=filename, iostat=io_error,position='append')
                write(1,*) length
                write(1,*)
                do i = 1, length
                    write(1,*) listofatoms(i)%name, &
                    listofatoms(i)%coords%x, &
                    listofatoms(i)%coords%y, &
                    listofatoms(i)%coords%z
                enddo
            close(1)
        end subroutine writetrjfile




        subroutine printmatrix(a)
            implicit none
            integer :: i,j
            real, dimension(:,:) :: a
            write(*,*)            
            do i = lbound(a,1), ubound(a,1)
               write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
            end do
        end subroutine printmatrix




        subroutine fileCleanUp(filename)
            implicit none
            logical :: file_exists
            character(len=72), intent(in) :: filename

            inquire(file=filename, exist=file_exists)
            if (file_exists .eqv. .true.) then
                open  (unit=1, file=filename, status="old")
                close (unit=1, status="delete")
            end if
        end subroutine fileCleanUp




        subroutine calcVelocities(listofFrames, nframes)!,velos)
            type( system ), dimension(:), allocatable, intent(inout) :: listofFrames
            integer, intent(in) :: nframes
            integer :: i,j,k
            do i = 1, nframes
                do j = 1, size(listofFrames(1)%subsystem) - 1
                    k = j + 1
                    listofFrames(i)%subsystem(j)%velo%v_x = ( listofFrames(i)%subsystem(k)%coords%x &
                                                            - listofFrames(i)%subsystem(j)%coords%x ) / 5.0d0
                    listofFrames(i)%subsystem(j)%velo%v_y = ( listofFrames(i)%subsystem(k)%coords%y &
                                                            - listofFrames(i)%subsystem(j)%coords%y ) / 5.0d0
                    listofFrames(i)%subsystem(j)%velo%v_z = ( listofFrames(i)%subsystem(k)%coords%z &
                                                            - listofFrames(i)%subsystem(j)%coords%z ) / 5.0d0

                    listofFrames(i)%subsystem(j)%s_velo = sqrt( (listofFrames(i)%subsystem(j)%velo%v_x)**2 &
                                                                + (listofFrames(i)%subsystem(j)%velo%v_y)**2 &
                                                                + (listofFrames(i)%subsystem(j)%velo%v_z)**2 )
                enddo
            enddo
        end subroutine calcVelocities
end module cp_mod