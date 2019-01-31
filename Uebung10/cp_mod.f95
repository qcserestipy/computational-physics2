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
            integer, dimension(:,:) :: a

            write(*,*)            
            do i = lbound(a,1), ubound(a,1)
               write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
            end do

        end subroutine printmatrix



        subroutine printrealmatrix(a)
            implicit none
            integer :: i,j
            real, dimension(:,:) :: a

            write(*,*)            
            do i = lbound(a,1), ubound(a,1)
               write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
            end do

        end subroutine printrealmatrix




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




        subroutine calcVelocities(listofFrames, nframes)
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




        real function distance(atomA, atomB)
            implicit none
            type( atom ) :: atomA, atomB

            distance = sqrt( (atomA%coords%x - atomB%coords%x)**2 &
                           + (atomA%coords%y - atomB%coords%y)**2 & 
                           + (atomA%coords%z - atomB%coords%z)**2 )
        end function distance




        subroutine ludcmps(A,L,U,n)
            implicit none
            real, dimension(:,:), allocatable, intent(in) :: A
            real, dimension(:,:), allocatable, intent(out) :: L
            real, dimension(:,:), allocatable, intent(out) :: U
            integer, intent(in) :: n
            real :: summ = 0.0d0
            integer :: i,j,k
            
            
            allocate(L(n,n))
            allocate(U(n,n))
            L = 0
            U = 0
            
            
            do j = 1, n  
              do i = j, n
                summ = 0.0d0
                do k = 1, j 
                  summ = summ + U(i,k) * L(k,j)
                enddo
                U(i,j) = A(i,j) - summ
              enddo
            
              do i = j, n 
                summ = 0.0d0
                do k = 1, j
                  summ = summ + U(j,k) * L(k,i)
                enddo
                if (U(j,j) == 0) then
                  write(*,*) "det(L) close to 0!\n Can't divide by 0..."
                  exit 
                end if
                L(j,i) = (A(j,i) - summ) / U(j,j)
              enddo
            enddo
        end subroutine ludcmps
            
            

            
        subroutine fwdsub(y_vec, l_mat, b_vec,n)
            implicit none
            real, intent(in), dimension(:,:), allocatable :: l_mat
            real, intent(out), dimension(:), allocatable :: y_vec
            real, intent(in), dimension(:), allocatable :: b_vec
            integer, intent(in) :: n
            real :: summ = 0.0d0
            integer :: i,j

            allocate(y_vec(n))

            y_vec=0
            y_vec(1) = b_vec(1) / l_mat(1,1)
            
            
            do i = 2, n
              y_vec(i) = b_vec(i) / l_mat(i,i)
              do j = 1, i - 1, 1
                summ = summ - l_mat(j,i) * y_vec(j)
              enddo
              summ = summ / l_mat(i,i)
              y_vec(i) = y_vec(i) + summ
              summ = 0.0d0
            enddo
            
        end subroutine fwdsub
            
            
        subroutine bwdsub(x_vec, u_mat, y_vec,n)
            implicit none
            real, intent(in), dimension(:,:), allocatable :: u_mat
            real, intent(in), dimension(:), allocatable :: y_vec
            real, intent(out),dimension(:), allocatable :: x_vec
            integer, intent(in) :: n
            integer :: i,j
            real :: summ = 0.0d0
            
            allocate(x_vec(n))
            x_vec(n) = y_vec(n) / u_mat(n,n)
            do i = n-1, 1, -1
              x_vec(i) = y_vec(i) / u_mat(i,i)
              do j = i + 1, n
                summ = summ - u_mat(j,i) * x_vec(j)
              enddo
              summ = summ / u_mat(i,i)
              x_vec(i) = x_vec(i) + summ
              summ = 0.0d0
            enddo

            write(*,*) "Solutionvector:"
            do i = 1, n
                write(*,*) x_vec(i)
            enddo
        end subroutine bwdsub




        subroutine tuple2int(i,j,k)
            implicit none
            integer, intent(in) :: i,j
            integer, intent(out) :: k
        
            if (i == 2 .and. j == 2) k =1
            if (i == 2 .and. j == 3) k =2
            if (i == 2 .and. j == 4) k =3
            if (i == 3 .and. j == 2) k =4
            if (i == 3 .and. j == 3) k =5
            if (i == 3 .and. j == 4) k =6
            if (i == 4 .and. j == 2) k =7
            if (i == 4 .and. j == 3) k =8
            if (i == 4 .and. j == 4) k =9
        end subroutine tuple2int
        
        
        
        
        subroutine int2tuple(k,i,j)
            implicit none
            integer, intent(out) :: i,j
            integer, intent(in) :: k
        
            if (k == 1) then
                i = 2 
                j = 2
            else if (k==2) then
                i = 2 
                j = 3
            else if (k==3) then
                i = 2 
                j = 4
            else if (k==4) then
                i = 3 
                j = 2
            else if (k==5) then
                i = 3 
                j = 3
            else if (k==6) then
                i = 3 
                j = 4
            else if (k==7) then
                i = 4 
                j = 2
            else if (k==8) then
                i = 4 
                j = 3
            else if (k==9) then
                i = 4 
                j = 4
            end if
        end subroutine int2tuple
        
        
        
        
        subroutine setUpGrid(pointMap,solutionVector)
            implicit none
            integer :: i,j, n, m, k
            integer, dimension(5,5), intent(out) :: pointMap
            real, dimension(:), allocatable, intent(out) :: solutionVector
            
            allocate(solutionVector(9))
            pointMap = 0
            solutionVector = 0
         
            do i = 0, 4
                pointMap(1,i+1) = i * 25 
            enddo
            do i = 1, 5
                pointMap(i+1,5) = pointMap(i,5) - 25
            enddo
        
            do i = 1, 9
                call int2tuple(i,n,m)
                solutionVector(i) = pointMap(n-1,m) + pointMap(n+1,m) + pointMap(n,m-1)+ pointMap(n,m+1)
            enddo
        
            do i = 2, 4
                do j = 2, 4
                    call tuple2int(i,j,k)
                    pointMap(i,j) = k
                enddo
            enddo 
        
        end subroutine setUpGrid
        
        
        
        
        subroutine setUpLgs(pointMap,matrix)
            implicit none 
            integer, dimension(5,5),intent(in) :: pointMap
            real, dimension(:,:), allocatable,intent(out) :: matrix
            integer :: i,j,k,l
        
            allocate(matrix(9,9))
            matrix = 0
            do i = 1, 9
                matrix(i,i) = 4.0d0
            enddo
        
            do i = 1,9
              do j = 1,9
                if(i == j + 1 .OR. i == j - 1 .OR. i == j + 3 .OR. i == j - 3) then
                   matrix(i,j) = -1.0d0
                end if
              enddo
            enddo            
        end subroutine
end module cp_mod