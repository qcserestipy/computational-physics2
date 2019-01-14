program gauss
IMPLICIT NONE
integer, parameter :: n=3
double precision a(n,n), b(n), x(n)
integer i,j
! matrix A
  data (a(1,i), i=1,3) /  4.0,  -9.0,  2.0 /
  data (a(2,i), i=1,3) /  2.0,  -4.0,  4.0 /
  data (a(3,i), i=1,3) /  -1.0,  2.0,  2.0 /
! matrix b
  data (b(i),   i=1,3) /  2.0,  3.0,  1.0 /

! print a header and the original equations
  write (*,200)
  do i=1,n
     write (*,201) (a(i,j),j=1,n), b(i)
  end do

  call gauss_jord(a,b,x,n)

! print matrix A and vector b after the elimination 
  write (*,202)
  do i = 1,n
     write (*,201)  (a(i,j),j=1,n), b(i)
  end do
! print solutions
  write (*,203)
  write (*,201) (x(i),i=1,n)
200 format (' Gauss elimination with scaling and pivoting ' &
    ,/,/,' Matrix A and vector b')
201 format (6f12.6)
202 format (/,' Matrix A and vector b after elimination')
203 format (/,' Solutions x(n)')
end


subroutine gauss_jord(a,b,x,n)
implicit none 
character(len=3) :: switch = 'aaa'
integer n
double precision a(n,n), b(n), x(n)
double precision s(n)
double precision c, pivot, store
integer i, j, k, l

! step 1: begin forward elimination
do k=1, n-1

! step 2: "scaling"
! s(i) will have the largest element from row i 
  do i=k,n                       ! loop over rows
    s(i) = 0.0
    do j=k,n                    ! loop over elements of row i
      s(i) = max(s(i),abs(a(i,j)))
    end do
  end do

write(*,*) "Should I use pivoting? (YES/no)"
read(*,*) switch
if (switch == 'yes' .or. switch == 'Yes' .or. switch == 'y' .or. switch == 'Y') then

write(*,*) "Pivoting!"
! step 3: "pivoting 1" 
! find a row with the largest pivoting element
    pivot = abs(a(k,k)/s(k))
    l = k
    do j=k+1,n
      if(abs(a(j,k)/s(j)) > pivot) then
        pivot = abs(a(j,k)/s(j))
        l = j
      end if
    end do

! Check if the system has a sigular matrix
    if(pivot == 0.0) then
      write(*,*) ' The matrix is sigular '
      return
    end if

! step 4: "pivoting 2" interchange rows k and l (if needed)
  if (l /= k) then
    do j=k,n
       store = a(k,j)
       a(k,j) = a(l,j)
       a(l,j) = store
    end do
    store = b(k)
    b(k) = b(l)
    b(l) = store
  end if
else 
  write(*,*) "Not pivoting!"

end if

write(*,*) "Starting elimination!"
! step 5: the elimination (after scaling and pivoting)
   do i=k+1,n
      c=a(i,k)/a(k,k)
      a(i,k) = 0.0
      b(i)=b(i)- c*b(k)
      do j=k+1,n
         a(i,j) = a(i,j)-c*a(k,j)
      end do
   end do
end do

! step 6: back substiturion 
x(n) = b(n)/a(n,n)
do i=n-1,1,-1
   c=0.0
   do j=i+1,n
     c= c + a(i,j)*x(j)
   end do 
   x(i) = (b(i)- c)/a(i,i)
end do

end subroutine gauss_jord

