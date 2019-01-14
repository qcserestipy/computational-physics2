!##################
!#      main      #
!##################
program euler
real, dimension(:), allocatable :: x_coord
real, dimension(:), allocatable :: y_coord
real :: h , ent_int
integer :: nsteps
character(len=72) :: output, output2

output = 'data001.dat'
output2 = 'dataI001.dat'

write(*,*) "Please type stepsize"
read(*,*) h

write(*,*) "Please type end of intervall"
read(*,*) ent_int

nsteps=int(ent_int/h)

allocate(x_coord(nsteps))
allocate(y_coord(nsteps))

call eulerSimple(x_coord, y_coord,h,nsteps)
call writeoutput(x_coord, y_coord,nsteps, output)

call eulerImproved(x_coord, y_coord,h,nsteps)
call writeoutput(x_coord, y_coord,nsteps, output2)

end program euler




!#################
!#   functions   #
!#################
subroutine eulerSimple(x_vec, y_vec, stepsize,nsteps)
real, dimension(nsteps), intent(out) :: x_vec
real, dimension(nsteps), intent(out) :: y_vec
integer, intent(in) :: nsteps
real, intent(in) :: stepsize
real :: y = 2.0d0, f = 0.0d0
real :: y_old, x = 0.0d0

do i = 1, nsteps
  x = x + stepsize
  y = y + stepsize * f
  f = f - stepsize * y_old
  x_vec(i) = x
  y_vec(i) = y
  y_old = y
enddo

end subroutine eulerSimple




subroutine eulerImproved(x_vec, y_vec, stepsize,nsteps)
real, dimension(nsteps), intent(out) :: x_vec
real, dimension(nsteps), intent(out) :: y_vec
integer, intent(in) :: nsteps
real, intent(in) :: stepsize
real :: y = 2.0d0, f = 0.0d0, x = 0.0d0, f_old

f_old = f
do i = 1, nsteps
  x = x + stepsize
  f = f - stepsize * y
  y = y + ( 0.5 * stepsize * ( f_old + f  )  )
  f_old = f
  x_vec(i) = x
  y_vec(i) = y
enddo
end subroutine eulerImproved




subroutine writeoutput(x_vec, y_vec, len, outfile)
real, dimension(len), intent(out) :: x_vec
real, dimension(len), intent(out) :: y_vec
integer :: len
character(len=*) :: outfile
open(unit=2,file=outfile)
  write(*,*) 'Schreibe Daten in outfile...'
  write(*,*) '-----------------------------------'
  do i = 1, len
    write(2,*) x_vec(i), y_vec(i)
  enddo
close(2)
end subroutine writeoutput
