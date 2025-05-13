! ============================================================================
! Name        : CP2.f90
! Author      : patrick
! Version     :
! Copyright   : Your copyright notice
! Description : Hello World in Fortran
! ============================================================================

program CP2
    integer :: imKreis = 0, ntests = 100000000
    real :: x,y, AbstandzumUrsprung, pi

    write(*,*) 'ntests', 'pi'
    write(*,*) '==============================='
!    do j = 1, 100000

        do i = 1, ntests
            call random_number(x)
            call random_number(y)
            AbstandzumUrsprung = SQRT(x**2 + y**2)
            if(AbstandzumUrsprung < 1.0) imKreis = imKreis + 1
        end do
        pi = 4.0*(REAL(imKreis)/REAL(ntests))
            write(*,*) ntests , pi


        ntests = ntests + 1000
        imKreis = 0

 !    end do
end program
