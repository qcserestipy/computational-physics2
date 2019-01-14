!########################
!#        main          #
!########################

program main 
    use cp_mod
    use veloverlet_mod


    !#########################
    !#       Aufgabe 2       #
    !#########################
    real :: dt = 0.02
    integer, parameter :: N = 120
    integer :: spec_points, i, j, io_error
    complex*16, dimension(:), allocatable :: f_x
    complex*16, dimension(:), allocatable :: f_k
    real*16, parameter :: PI = 4 * ATAN(1.)
    real, dimension(:), allocatable :: t, power_spec, dummy, omega
    character(len=72) :: infile = 'CO.xyz'
    type ( atom ), dimension(:), allocatable :: atoms

    call readxyzfile(infile, atoms)
    call timeloop(atoms, N, dt)

    ! Allocate x-space data points
    spec_points = int(N/dt) 
    allocate(t(spec_points))
    allocate(omega(spec_points))
    allocate(dummy(spec_points))
    allocate(power_spec(spec_points))
    allocate(f_k(spec_points))
    allocate(f_x(spec_points))

    open(unit=1,file='velocity.out')
        do i = 1, spec_points
            read(1,*) t(i), dummy(i)
            f_x(i) = cmplx(0,dummy(i))
        enddo
    close(1)        

   
  
    ! Execute the discrete fourier transform
    ! It generates N data points in k-space,
    ! and they are assigned to the array f_k
    call dft(f_x, f_k, spec_points)
    
    do i = 1, spec_points
        omega(i) = i * ((2*PI)/spec_points*dt)
    enddo

    write(*,*) "####################################"
    write(*,*) "#        Entering power spec       #" 
    write(*,*) "####################################"
    open(unit=1,file="power_spec.out")
    do i = 1, spec_points
        power_spec(i) = abs(f_k(i))**2
        write(1,*) omega(i), power_spec(i)
    enddo
    close(1)

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
    implicit none
    
    integer, intent(in) :: N
    complex*16, dimension(N), intent(in) :: f_x
    complex*16, dimension(N), intent(out) :: f_k
    real*16, parameter :: PI = 4 * ATAN(1.)
    integer :: i,j

    write(*,*) "####################################"
    write(*,*) "#       Entering F-transform       #" 
    write(*,*) "####################################"
    open(unit=1,file="FT.out")
    do i=1,N
      ! initializing f_k(i)
      f_k(i) = 0
      
      ! main calculation part
      
      do j=1,N
          f_k(i) = f_k(i) + f_x(j) * cexp(cmplx(0, -2*PI*i*j/N))
          
      end do
      write(1,*) i, f_k(i)
    end do
    close(1)

end subroutine dft


! Super simple routine that executes discrete fourier back transform
! f_x : x-space data points (double precision complex array with size N)
! f_k : k-space data points (double precision complex array with size N)
! N   : the size of f_x and f_k
subroutine bdft(f_x, f_k, N)
    implicit none
    integer, intent(in) :: N
    complex*16, dimension(N), intent(out) :: f_x
    complex*16, dimension(N), intent(in) :: f_k
    real*16, parameter :: PI = 4 * ATAN(1.)
    integer :: i,j

    do i=1,N
      ! initializing f_k(i)
      f_x(i) = 0
      
      ! main calculation part
      do j=1,N
          f_x(i) = f_k(i) + f_k(j) * cexp(cmplx(0, 2*PI*i*j/N))
      end do
    end do

end subroutine bdft