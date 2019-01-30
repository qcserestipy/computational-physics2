program LU
real, dimension(3,3) :: L_mat(3,3)
real, dimension(3,3) :: U_mat(3,3)
real, dimension(3,3) :: A_mat(3,3)
real, dimension(3) :: b_vec(3)
real, dimension(3) :: y_vec(3)
real, dimension(3) :: x_vec(3)
! matrix A
 A_mat = reshape( (/4.0,  -9.0,  2.0, 2.0,  -4.0,  4.0, -1.0,  2.0,  2.0/), (/3, 3/) )
! matrix b
 b_vec = (/2.0,  3.0,  1.0/)


write(*,*) "===================="
write(*,*) "Matrix A"
write(*,*) "===================="
write(*,*) A_mat(1,1), A_mat(2,1), A_mat(3,1)
write(*,*) A_mat(1,2), A_mat(2,2), A_mat(3,2)
write(*,*) A_mat(1,3), A_mat(2,3), A_mat(3,3)
write(*,*) "===================="


write(*,*) "Vector b"
write(*,*) "===================="
write(*,*) b_vec(1)
write(*,*) b_vec(2)
write(*,*) b_vec(3)
write(*,*) "===================="



CALL ludcmps(A_mat,L_mat,U_mat,3)

write(*,*) "Matrix L"
write(*,*) "===================="
write(*,*) L_mat(1,1), L_mat(2,1), L_mat(3,1)
write(*,*) L_mat(1,2), L_mat(2,2), L_mat(3,2)
write(*,*) L_mat(1,3), L_mat(2,3), L_mat(3,3)
write(*,*) "===================="

write(*,*) "Matrix U"
write(*,*) "===================="
write(*,*) U_mat(1,1), U_mat(2,1), U_mat(3,1)
write(*,*) U_mat(1,2), U_mat(2,2), U_mat(3,2)
write(*,*) U_mat(1,3), U_mat(2,3), U_mat(3,3)
write(*,*) "===================="

CALL fwdsub(y_vec, L_mat, b_vec,3)
write(*,*) "Vector y"
write(*,*) "===================="
write(*,*) y_vec(1)
write(*,*) y_vec(2)
write(*,*) y_vec(3)
write(*,*) "===================="

CALL bwdsub(x_vec, U_mat, y_vec,3)
write(*,*) "Vector x"
write(*,*) "===================="
write(*,*) x_vec(1)
write(*,*) x_vec(2)
write(*,*) x_vec(3)
write(*,*) "===================="


write(*,*) "approximate determinant of matrix A: ", det_approx(U_mat,3)
write(*,*) "exact determinant of matrix A: ", 8


end program LU


subroutine ludcmps(A,L,U,n)
real, dimension(:,:), allocatable, intent(in) :: A
real, dimension(:,:), allocatable, intent(out) :: L
real, dimension(:,:), allocatable, intent(out) :: U
real :: summ = 0.0d0

do i = 1, n
  do j = 1,n 
    L(i,j) = 0.0d0
    U(i,j) = 0.0d0
  enddo
enddo


do j = 1, n  
  do i = j, n
    summ = 0.0d0
    do k = 1, j 
      summ = summ + U(i,k) * L(k,j)
    enddo
    U(i,j) = A(i,j) - summ
  enddo

  do i = j, n 
    summ = 0.0d0
    do k = 1, j
      summ = summ + U(j,k) * L(k,i)
    enddo
    if (U(j,j) == 0) then
      write(*,*) "det(L) close to 0!\n Can't divide by 0..."
      exit 
    end if
    L(j,i) = (A(j,i) - summ) / U(j,j)
  enddo
enddo


end subroutine ludcmps



subroutine fwdsub(y_vec, l_mat, b_vec,n)
real, intent(in), dimension(:,:), allocatable :: l_mat
real, intent(out),dimension(:), allocatable :: y_vec
real, intent(in),dimension(:), allocatable :: b_vec
real :: summ = 0.0d0

allocate(u_mat(n,n))
allocate(y_vec(n))
allocate(b_vec(n))
y_vec=0
y_vec(1) = b_vec(1) / l_mat(1,1)


do i = 2, n
  y_vec(i) = b_vec(i) / l_mat(i,i)
  do j = 1, i - 1, 1
    summ = summ - l_mat(j,i) * y_vec(j)
  enddo
  summ = summ / l_mat(i,i)
  y_vec(i) = y_vec(i) + summ
  summ = 0.0d0
enddo

end subroutine fwdsub


subroutine bwdsub(x_vec, u_mat, y_vec,n)
  real, intent(in), dimension(:,:), allocatable :: u_mat
  real, intent(out),dimension(:), allocatable :: y_vec
  real, intent(in),dimension(:), allocatable :: b_vec
real :: summ = 0.0d0

allocate(u_mat(n,n))
allocate(y_vec(n))
allocate(b_vec(n))

x_vec(n) = y_vec(n) / u_mat(n,n)
do i = n-1, 1, -1
  x_vec(i) = y_vec(i) / u_mat(i,i)
  do j = i + 1, n
    summ = summ - u_mat(j,i) * x_vec(j)
  enddo
  summ = summ / u_mat(i,i)
  x_vec(i) = x_vec(i) + summ
  summ = 0.0d0
enddo
end subroutine bwdsub

function det_approx(array,n)
real, intent(in), dimension(3,3) :: array(3,3)
real :: temp = 1.0d0
do i = 1, n
  temp = temp * array(i,i)
enddo 

det_approx = temp
return
end function det_approx


