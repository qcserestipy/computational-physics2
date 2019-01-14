program euler_method
use iso_fortran_env, only: real64
implicit none
 
abstract interface
  ! a derivative dy/dt as function of y and t
  function derivative(y, t)
    use iso_fortran_env, only: real64
    real(real64) :: derivative
    real(real64), intent(in) :: t, y
  end function
end interface
 
real(real64), parameter :: a = 0, b = 100, &
    h(4) = [0.001, 0.01, 0.02, 0.05]
 
integer :: i

call euler(0.0d0,2.0d0,0,100,h(1)) 
 
contains
 
subroutine euler(f, y0, a, b, h)
  procedure(derivative) :: f
  real(real64), intent(in) :: y0, a, b, h
  real(real64) :: t, y
 
  if (a > b) return
  if (h <= 0) stop "negative step size"
 
  print '("# h = ", F0.3)', h
 
  y = y0
  t = a
 
  do
    print *, t, y
    t = t + h
    if (t > b) return
    y = y + h * f(y, t)
  end do
end subroutine

end program 
