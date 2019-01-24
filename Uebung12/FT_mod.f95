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
        
    subroutine autoCorrelationCalculator(list,corrA,corrB,corrC,corrG)
        implicit none
        real, dimension(18000), intent(in) :: list 
        real, dimension(6000) :: vx
        real, dimension(6000) :: vy
        real, dimension(6000) :: vz
        real, dimension(2000) :: vA
        real, dimension(2000) :: vB
        real, dimension(2000) :: vC
        real, dimension(2000) :: vG
        real, dimension(2000), intent(out) :: corrA
        real, dimension(2000), intent(out) :: corrB
        real, dimension(2000), intent(out) :: corrC
        real, dimension(2000), intent(out) :: corrG
        integer :: i
        
        do i = 1, 18000, 3
            vx(i/3) = list(i)
        enddo

        do i = 2, 18000, 3
            vy(i/3) = list(i)
        enddo

        do i = 3, 18000, 3
            vz(i/3) = list(i)
        enddo

        do i = 1, 6000, 3
            vA(i/3) = sqrt( vx(i)**2 + vy(i)**2 + vz(i)**2 )
        enddo

        do i = 2, 6000, 3
            vB(i/3) = sqrt( vx(i)**2 + vy(i)**2 + vz(i)**2 )
        enddo

        do i = 3, 6000, 3
            vC(i/3) = sqrt( vx(i)**2 + vy(i)**2 + vz(i)**2 )
        enddo

        do i = 1, 2000
            vG(i) = (vA(i) + vB(i) + vC(i))/3
        enddo

        open(unit=1,file="autocorr.out")
        do i = 1, 1999
            corrA(i) =  vA(i)*vA(i+1)
            corrB(i) =  vB(i)*vB(i+1)
            corrC(i) =  vC(i)*vC(i+1)
            corrG(i) =  vG(i)*vG(i+1)
            write(1,*) i, corrA(i), corrB(i), corrC(i), corrG(i)
        enddo
        close(1)


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