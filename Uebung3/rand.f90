program lcg

real*4 ::  zeta=0, oldzeta=0, first_moment=0, second_moment=0, corr=0
integer :: c=5701, a=3612, m=566927, cycler = 0
character(len=1), parameter :: TAB = achar(9) 

do n = 1, 10000
  call random_number(zeta)
  corr = corr + oldzeta*zeta
  oldzeta = zeta
  ! write output
  open (unit=20,file="plot1.dat",position='append')
    write(20,*) n, zeta
  close(20)

  write(*,*) n, zeta
  first_moment = first_moment + zeta
  second_moment = second_moment + zeta**2
  
  cycler = cycler + 1
enddo

  first_moment = first_moment / REAL(cycler)
  second_moment = second_moment / REAL(cycler)
  corr = corr / REAL(cycler)
  write(*,*) "first moment", TAB , "second moment",TAB, "correlation"
  write(*,*) first_moment, second_moment, corr
end program lcg


