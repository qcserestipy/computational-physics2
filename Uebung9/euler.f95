!##################
!#      main      #
!##################
program euler
real, dimension(:), allocatable :: x_coord
real, dimension(:), allocatable :: y_coord
real :: h , ent_int
integer :: nsteps


write(*,*) "Please type stepsize"
read(*,*) h

write(*,*) "Please type end of intervall"
read(*,*) ent_int

nsteps=int(ent_int/h)

allocate(x_coord(nsteps))
allocate(y_coord(nsteps))


call eulerImproved(x_coord, y_coord,h,nsteps)


end program euler




!#################
!#   functions   #
!#################


subroutine eulerImproved(x_vec, y_vec, stepsize,nsteps)
real, dimension(nsteps), intent(out) :: x_vec
real, dimension(nsteps), intent(out) :: y_vec
integer, intent(in) :: nsteps
real, intent(in) :: stepsize
real :: y = 2.3d0, f = 0.0d0, x = 0.0d0, f_old, E = 0.0d0

f_old = f
open(unit=2,file='CO.xyz')
open(unit=3,file='E.out')
open(unit=4,file='period.out')
write(*,*) '====================='
write(*,*) '    nsteps    ', 'distance          ' , 'velocity    ', 'force'
do i = 1, nsteps
  write(2,*) '2'
  write(2,*) 
  write(2,*) 'C', 0.0, 0.0, 0.0
  write(2,*) 'O', 0.0, 0.0, abs(y)
  E = updateEnergy(abs(y))
  write(3,*) abs(y) , E
  x = x + stepsize
  f = f - stepsize * y
  y = y + ( 0.5 * stepsize * ( f_old + f  )  )
  f_old = f
  x_vec(i) = x
  y_vec(i) = y
enddo
close(4)
close(3)
close(2)
end subroutine eulerImproved

function updateEnergy(x)
  real :: x
  real, parameter :: D = 0.493172, a = 2.05666, R_0 = 1.20106
      updateEnergy = D*(( 1 - exp(-a*(x - R_0)) )**2)
  end function updateEnergy
