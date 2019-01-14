!##################
!#      main      #
!##################
program euler
!real, dimension(:), allocatable :: x_coord
!real, dimension(:), allocatable :: y_coord
!real, dimension(:), allocatable :: t
real :: h , ent_int
integer :: nsteps
character(len=72) :: output

output = 'planets.xyz'

write(*,*) "Please type stepsize"
read(*,*) h

write(*,*) "Please type end of intervall"
read(*,*) ent_int

nsteps=int(ent_int/h)

!allocate(x_coord(nsteps))
!allocate(y_coord(nsteps))
!allocate(t(nsteps))

call eulerImproved(h,nsteps)
!call writeoutput(nsteps, output)

end program euler




!#################
!#   functions   #
!#################

subroutine eulerImproved(stepsize,nsteps)
!real, dimension(nsteps), intent(out) :: x_vec
!real, dimension(nsteps), intent(out) :: y_vec
!real, dimension(nsteps), intent(out) :: t_vec
integer, intent(in) :: nsteps
real, intent(in) :: stepsize
real :: y = 0.0d0, f_x = 0.0d0, f_y = 1.63d0, x = 0.5d0, f_xold, f_yold, r = 1.0d0, x_old, y_old, r_old

open(unit=2,file='lol.xyz')
  write(*,*) 'Schreibe Daten in outfile...'
  write(*,*) '-----------------------------------'
  
x_old=x
y_old=y
f_xold = f_x
f_yold = f_y
do i = 1, nsteps
  t = t + stepsize
  
  r_old = SQRT(x_old**2 + y_old**2)
  
  !Prediktor-Schritt
  x = x_old + stepsize * f_xold
  y = y_old + stepsize * f_yold
  f_y = f_yold - (stepsize * y_old)*(1/(r_old**3.5))
  f_x = f_xold - (stepsize * x_old)*(1/(r_old**3.5))

  r = SQRT(x**2 + y**2)
  
  !Korrektor-Schritt
  y = y_old + ( 0.5 * stepsize * ( f_yold + f_y  )  )
  x = x_old + ( 0.5 * stepsize * ( f_xold + f_x  )  )
  f_y = f_yold + ( 0.5 * stepsize * ( -y_old/(r_old**3.5) - y/(r**3.5)  )  )
  f_x = f_xold + ( 0.5 * stepsize * ( -x_old/(r_old**3.5) - x/(r**3.5)  )  )
  
  
  f_yold = f_y
  f_xold = f_x
  x_old=x
  y_old=y
  !x_vec(i) = x
  !y_vec(i) = y
  
  !t_vec(i) = t
  write(2,*) '2'
  write(2,*)
  write(2,*) 'H', x, y, 0.0
  write(2,*) 'Au', 0.0, 0.0, 0.0

enddo
close(2)
end subroutine eulerImproved




!subroutine writeoutput(x_vec, y_vec, t_vec, len, outfile)
!real, dimension(len), intent(out) :: x_vec
!real, dimension(len), intent(out) :: y_vec
!real, dimension(len), intent(out) :: t_vec
!integer :: len
!character(len=*) :: outfile
!!open(unit=2,file=outfile)
!  write(*,*) 'Schreibe Daten in outfile...'
!  write(*,*) '-----------------------------------'
!  write(2,*) '2'
!  write(2,*)
!  do i = 1, len
!    write(2,*) 'H', x_vec(i), y_vec(i), 0.0
!    write(2,*) 'Au', 0.0, 0.0, 0.0
!  enddo
!close(2)
!end subroutine writeoutput
