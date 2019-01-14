! ============================================================================
! Name        : importance_sampling.f90
! Author      : patrick
! Version     :
! Copyright   : Your copyright notice
! Description : importance sampling program
! ============================================================================

program mWMC
    real :: f, w, x, u, I = 0.0
    integer :: ntests = 10000


    do j = 1, ntests
        call random_number(u)
        x = xfunc(u)
        w = wfunc(x)
        f = func(x)
        I = I + (f/w)
    end do
    I = I / REAL(ntests)
    write(*,*) ntests, I


end program mWMC

function func(i) result(j) !berechnet die Funktion f(x)
    real, intent(in) :: i
    real             :: j
    j = i**(-1.0/3.0) + (i/10.0)
end function func

function xfunc(i) result(j) !berechnet die Inverse von u(x)
    real, intent(in) :: i
    real             :: j
    j = ( (i+0.2746) / 1.2746)**(3.0/2.0)
end function xfunc

function wfunc(i) result(j) !berechnet die Gewichtungsfunktion w(x)
    real, intent(in) :: i
    real             :: j
    j = 0.8497*i**(-1.0/3.0)
end function wfunc

