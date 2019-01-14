!###############
!#    main     #
!###############
program veloverlet
real :: x = 2.3
real :: v = 0, dt = 0.02d0, f = 0.0d0, m = 15.9994, E = 0.0d0, t = 0
integer :: nsteps
nsteps = int(1000/dt)
!###################
!#    timeloop     #
!###################

open(unit=2,file='CO.xyz')
open(unit=3,file='E.out')
open(unit=4,file='period.out')
write(*,*) '====================='
write(*,*) '    nsteps    ', 'distance          ' , 'velocity    ', 'force'
do i = 1, nsteps
  write(2,*) '2'
  write(2,*) 
  write(2,*) 'C', 0.0, 0.0, 0.0
  write(2,*) 'O', 0.0, 0.0, x
  write(3,*) x, E
  write(4,*) t, x
  v = v + f * dt / (2*m)
  x = x + v * dt
  f =  updateForce(x)
  E = updateEnergy(x)
  t = t + dt
  write(*,*) i , x, v, f
enddo
close(4)
close(3)
close(2)
end program veloverlet

!###################
!#    functions    #
!###################

function updateForce(x)
real :: x
real, parameter :: D = 0.493172, a = 2.05666, R_0 = 1.20106
    updateForce = - 2 * a * D * exp(a*(R_0 - x)) * (1 - exp(a*(R_0 - x)))
end function updateForce


function updateEnergy(x)
real :: x
real, parameter :: D = 0.493172, a = 2.05666, R_0 = 1.20106
    updateEnergy = D*(( 1 - exp(-a*(x - R_0)) )**2)
end function updateEnergy
