! ============================================================================
! Name        : mittelwert-monte-carlo.f90
! Author      : patrick
! Version     :
! Copyright   : Your copyright notice
! Description : Hello World in Fortran
! ============================================================================

program mWMC
    real :: f = 0.0, x, b = 1.0, a = 0.1, I, I_anal
    integer :: nhits, ntests = 100

    write(*,*) 'ntests     ', 'num. Int'
    do j = 1, 6

    do i = 1, ntests
        call random_number(x)
        f = f + x**(-1.0/3.0) + (x/10.0)
    end do

    I = ((b-a)/ntests) * f

    write(*,*) ntests, I
    ntests = ntests * 10
    nhits = 0
    end do

    I_anal = func(b) - func(a)
    write(*,*) I_anal

end program mWMC

function func(i) result(j) !berechnet die Stammfunktion
    real, intent(in) :: i
    real             :: j
    j = (3.0/2.0)*i**(2.0/3.0) + ((i**2) / 20.0)
end function func

