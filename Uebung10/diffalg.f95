program main
    use cp_mod

    implicit none
    integer, dimension(5,5) :: pointMap
    integer, dimension(9) :: solutionVector
    integer, dimension(9,9) :: matrix
    integer :: i,j
   
    call setUpGrid(pointMap,solutionVector)
    call printmatrix(pointMap)
    call setUpLgs(pointMap,matrix)
    !call printmatrix(matrix)
    
    CALL ludcmps(matrix,L_mat,U_mat,9)

    CALL fwdsub(y_vec, L_mat, b_vec,9)

    CALL bwdsub(x_vec, U_mat, y_vec,9)
end program main




subroutine tuple2int(i,j,k)
    implicit none
    integer, intent(in) :: i,j
    integer, intent(out) :: k

    if (i == 2 .and. j == 2) k =1
    if (i == 2 .and. j == 3) k =2
    if (i == 2 .and. j == 4) k =3
    if (i == 3 .and. j == 2) k =4
    if (i == 3 .and. j == 3) k =5
    if (i == 3 .and. j == 4) k =6
    if (i == 4 .and. j == 2) k =7
    if (i == 4 .and. j == 3) k =8
    if (i == 4 .and. j == 4) k =9
end subroutine tuple2int




subroutine int2tuple(k,i,j)
    implicit none
    integer, intent(out) :: i,j
    integer, intent(in) :: k

    if (k == 1) then
        i = 2 
        j = 2
    else if (k==2) then
        i = 2 
        j = 3
    else if (k==3) then
        i = 2 
        j = 4
    else if (k==4) then
        i = 3 
        j = 2
    else if (k==5) then
        i = 3 
        j = 3
    else if (k==6) then
        i = 3 
        j = 4
    else if (k==7) then
        i = 4 
        j = 2
    else if (k==8) then
        i = 4 
        j = 3
    else if (k==9) then
        i = 4 
        j = 4
    end if
end subroutine int2tuple




subroutine setUpGrid(pointMap,solutionVector)
    implicit none
    integer :: i,j, n, m, k
    integer, dimension(5,5), intent(out) :: pointMap
    integer, dimension(9), intent(out) :: solutionVector
    
    pointMap = 0
    solutionVector = 0
 
    do i = 0, 4
        pointMap(1,i+1) = i * 25 
    enddo
    do i = 1, 5
        pointMap(i+1,5) = pointMap(i,5) - 25
    enddo

    do i = 1, 9
        call int2tuple(i,n,m)
        solutionVector(i) = pointMap(n-1,m) + pointMap(n+1,m) + pointMap(n,m-1)+ pointMap(n,m+1)
    enddo

    do i = 2, 4
        do j = 2, 4
            call tuple2int(i,j,k)
            pointMap(i,j) = k
        enddo
    enddo 

end subroutine setUpGrid




subroutine setUpLgs(pointMap,matrix)
    use cp_mod
    implicit none 
    integer, dimension(5,5) :: pointMap
    integer, dimension(9,9) :: matrix
    integer :: i,j,k,l

    matrix = 0
    do i = 1, 9
        matrix(i,i) = 4
    enddo

    do i = 1,9
      do j = 1,9
        if(i == j + 1 .OR. i == j - 1 .OR. i == j + 3 .OR. i == j - 3) then
           matrix(i,j) = -1
        end if
      enddo
    enddo

    call printmatrix(matrix)
    
end subroutine