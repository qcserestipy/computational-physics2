! als x Intervall wurde [0.1,1] verwendet weshalb die zufallszahl für x
!neugewürfelt wird wenn der wert kleiner ist als x
!für das y intervall [0,2.17] wurde das globale maximum bei 0.1 bestimmt
!die Fläche berechent sich dann aus 2.17*0.9=1.953
program nm_int
    real :: x, y, I, area = 1.953
    integer :: nhits=0, ntests = 100

 write(*,*) 'ntests     ', 'num. Int'
    do j = 1, 6
    do i = 1, ntests
        do while(x < 0.1)
            call random_number(x)
        end do
            call random_number(y)
            y = y * 2.17

       if (y < func(x)) nhits = nhits + 1
    end do

    I = (REAL(nhits)/REAL(ntests)) * area
    write(*,*) ntests, I
    ntests = ntests * 10
    nhits = 0
    end do




end program nm_int

function func(i) result(j) !berechnet den Funktionswert der Stelle i
    real, intent(in) :: i
    real             :: j
    j = i**(-1.0/3.0) + i/10.0
end function func
