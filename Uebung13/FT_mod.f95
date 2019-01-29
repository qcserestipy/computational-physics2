module FT_mod
    use cp_mod
    implicit none
    save


    !##############################
    !#    function definitions    #
    !##############################
    contains

    subroutine MDReader(filename, t_velos)
        implicit none
        character(len=72),intent(in) :: filename
        real, dimension(18000), intent(inout) :: t_velos
        integer :: i
        real :: dummy
        open(unit=1,file=filename)
        do i = 1, 18000, 3
            read(1,*) dummy, dummy, dummy, dummy, t_velos(i), t_velos(i+1), t_velos(i+2)
        enddo
        close(1)
    end subroutine MDReader
    
    


    subroutine autoCorrelationCalculator(listofFrames, nframes,correlation)
        implicit none
        type( system ), dimension(:), allocatable, intent(in) :: listofFrames
        real, dimension(:,:), allocatable, intent(out) :: correlation
        integer, intent(in) :: nframes
        integer :: i,j,k, natoms

        natoms = size(listofFrames(1)%subsystem)
        allocate(correlation(nframes,natoms))

        do i = 1, nframes - 1
            do j = 1, natoms
                k = i + 1
                correlation(i,j) = (real(listofFrames(i)%subsystem(j)%s_velo) &
                                  * real(listofFrames(k)%subsystem(j)%s_velo) )
            enddo
        enddo
    end subroutine autoCorrelationCalculator




    subroutine discreteFourierTransformer(f_x, f_k, N)
        implicit none
        integer, intent(in) :: N
        complex*16, dimension(N), intent(in) :: f_x
        complex*16, dimension(N), intent(out) :: f_k
        real*16, parameter :: PI = 4 * ATAN(1.)
        integer :: i,j

        write(*,*) "####################################"
        write(*,*) "#       Entering F-transform       #" 
        write(*,*) "####################################"

        do i=1,N
            f_k(i) = 0
            do j=1,N
                f_k(i) = f_k(i) + f_x(j) * cexp(cmplx(0, -2*PI*i*j/N))
            end do
            write(*,*) i, f_k(i)
        end do
    end subroutine discreteFourierTransformer




    subroutine powerSpectrum(dataPoints,spec_points,outfile)
        implicit none
        character(len=72), intent(in) :: outfile
        integer, intent(in) :: spec_points
        real, dimension(spec_points), intent(in) :: dataPoints 
        complex*16, dimension(:), allocatable :: f_x
        complex*16, dimension(:), allocatable :: f_k
        real*16, parameter :: PI = 4 * ATAN(1.)
        real, dimension(:), allocatable :: power_spec, omega
        real :: dt = 0.002
        integer :: i

        allocate(omega(spec_points))
        allocate(power_spec(spec_points))
        allocate(f_k(spec_points))
        allocate(f_x(spec_points))

        do i = 1, spec_points
            f_x(i) = cmplx(0,dataPoints(i))
        enddo

        call discreteFourierTransformer(f_x, f_k, spec_points)

        do i = 1, spec_points
            omega(i) = i * ((2*PI)/(spec_points*dt))
        enddo

        write(*,*) "#########################################"
        write(*,*) "#        Entering power spectrum        #" 
        write(*,*) "#########################################"
        open(unit=1,file=outfile)
        do i = 1, spec_points
            power_spec(i) = abs(f_k(i))**2
            write(1,*) omega(i), power_spec(i)
        enddo
        close(1)
    end subroutine powerSpectrum


end module FT_mod