!########################
!#        main          #
!########################

program main 
    use cp_mod
    use veloverlet_mod


    !#########################
    !#       Aufgabe 1       #
    !#########################
    real :: t_0 = -1.2, t_N = 1.2, dt = 0.02
    integer, parameter :: N = 120
    complex*16, dimension(N) :: f_x
    complex*16, dimension(N) :: f_k
    integer i,j
    character(len=72) :: infile = 'CO.xyz'
    character(len=72) :: outfile = 'output.txt'
    type ( atom ), dimension(:), allocatable :: atoms

    call readxyzfile(infile, atoms)
    OPEN  (UNIT=1, FILE="CO.trj", STATUS="OLD")
    CLOSE (UNIT=1, STATUS="DELETE")
    call timeloop(atoms, 100, 0.02)


    ! Allocate x-space data points
    do i=1,N
      f_x(i) = y_k(t_0)
      t_0 = t_0 + dt
    end do
    ! Execute the discrete fourier transform
    ! It generates N data points in k-space,
    ! and they are assigned to the array f_k
    write(*,*) "######################### transformation #########################"
    call dft(f_x, f_k, N)

    ! Check
    do i=1,N
      write(*,*) f_k(i) 
    end do

write(*,*) "######################### backtransformation #########################"
    !back transformation
    call bdft(f_x, f_k, N)

    ! Check
    do i=1,N
        write(*,*) f_x(i) 
      end do
end program main

!########################
!#      functions       #
!########################


real function y_k(x)
    real :: x
    
    if (abs(x) < 1) then
        y_k = 1
    else if (abs(x) >= 1) then
        y_k = 0
    end if

end function y_k


! Super simple routine that executes discrete fourier transform
! f_x : x-space data points (double precision complex array with size N)
! f_k : k-space data points (double precision complex array with size N)
! N   : the size of f_x and f_k
subroutine dft(f_x, f_k, N)
    integer, intent(in) :: N
    complex*16, dimension(N), intent(in) :: f_x
    complex*16, dimension(N), intent(out) :: f_k
    real*16, parameter :: PI = 4 * ATAN(1.)

    integer i,j
    do i=1,N
      ! initializing f_k(i)
      f_k(i) = 0
      
      ! main calculation part
      do j=1,N
          f_k(i) = f_k(i) + f_x(j) * cexp(cmplx(0, -2*PI*i*j/N))
      end do
    end do

end subroutine dft


! Super simple routine that executes discrete fourier back transform
! f_x : x-space data points (double precision complex array with size N)
! f_k : k-space data points (double precision complex array with size N)
! N   : the size of f_x and f_k
subroutine bdft(f_x, f_k, N)
    integer, intent(in) :: N
    complex*16, dimension(N), intent(out) :: f_x
    complex*16, dimension(N), intent(in) :: f_k
    real*16, parameter :: PI = 4 * ATAN(1.)

    integer i,j
    do i=1,N
      ! initializing f_k(i)
      f_x(i) = 0
      
      ! main calculation part
      do j=1,N
          f_x(i) = f_k(i) + f_k(j) * cexp(cmplx(0, 2*PI*i*j/N))
      end do
    end do

end subroutine bdft