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
        type( velocity ), dimension(6000), intent(inout) :: t_velos
        integer :: i
        real :: dummy

        open(unit=1,file=filename)
        do i = 1, 6000
            read(1,*) dummy, dummy, dummy, dummy, t_velos(i)%v_x, t_velos(i)%v_y, t_velos(i)%v_z
        enddo
        close(1)
    end subroutine MDReader

    subroutine autoCorrelation(outfile,velos)
        implicit none
        real, dimension(6000) :: correlation

    end subroutine autoCorrelation


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


end module FT_mod