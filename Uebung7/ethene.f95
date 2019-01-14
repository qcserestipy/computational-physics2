Program Main
implicit none
integer, parameter :: n=2
double precision :: a(n,n), x(n,n)
double precision, parameter:: abserr=1.0e-09
integer i, j

! matrix A
  data (a(1,i), i=1,2) /   0.0,  1.0 /
  data (a(2,i), i=1,2) /   1.0,  0.0 /
  
! print a header and the original matrix
  write (*,200)
  do i=1,n
     write (*,201) (a(i,j),j=1,n)
  end do

  call Jacobi(a,x,abserr,n)

! print solutions
  write (*,202)
  write (*,201) (a(i,i),i=1,n)
  write (*,203)
  do i = 1,n
     write (*,201)  (x(i,j),j=1,n)
  end do

200 format (' Eigenvalues and eigenvectors (Jacobi method) ',/, &
            ' Matrix A')
201 format (6f12.6)
202 format (/,' Eigenvalues')
203 format (/,' Eigenvectors')
end program main







subroutine Jacobi(a,x,abserr,n)
implicit none
integer i, j, k, n
double precision a(n,n),x(n,n)
double precision abserr, b2, bar
real :: summ = 0.d0
double precision beta, coeff, c, s, cs, sc
integer :: nsweeps = 0
! initialize x(i,j)=0, x(i,i)=1
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
      if (i .ne. j) summ = summ + abs(a(i,j))**2
      nsweeps = nsweeps + 1
    end do
  end do
end do
write(*,*) "=========================================================="
write(*,*) "Sum of squared non diagonal elements and number of sweeps"
write(*,*) nsweeps, summ
write(*,*) "=========================================================="

return
end
