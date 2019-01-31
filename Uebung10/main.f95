program main
    use cp_mod

    implicit none
    integer, dimension(5,5) :: pointMap
    real, dimension(5,5) :: heatMap
    real, dimension(:), allocatable :: solutionVector
    real, dimension(:,:), allocatable :: matrix
    integer :: i,j,k
    real, dimension(:,:), allocatable :: L_mat
    real, dimension(:,:), allocatable :: U_mat
    real, dimension(:), allocatable :: y_vec
    real, dimension(:), allocatable :: x_vec
   
    write(*,*) "####################################"
    write(*,*) "#     Setting up grid on plate     #" 
    write(*,*) "####################################"
    call setUpGrid(pointMap,solutionVector)
    call printmatrix(pointMap)


    write(*,*) "####################################"
    write(*,*) "#         Setting up LGS           #" 
    write(*,*) "####################################"
    call setUpLgs(pointMap,matrix)
    call printrealmatrix(matrix)

    write(*,*) "####################################"
    write(*,*) "#    Entering LU decomposition     #" 
    write(*,*) "####################################"
    call ludcmps(matrix,L_mat,U_mat,9)
    call fwdsub(y_vec, L_mat, solutionVector,9)
    call bwdsub(x_vec, U_mat, y_vec,9)

    write(*,*) "####################################"
    write(*,*) "#    Displaying final heat map     #" 
    write(*,*) "####################################"

    heatMap = real(pointMap)
    do i = 2, 4
        do j = 2, 4
            call tuple2int(i,j,k)
            heatMap(i,j) = x_vec(k)
        enddo
    enddo
    call printrealmatrix(heatMap)

    open(unit=1,file="heatMap.out")
    do j = 1, 5
        write(1,*) heatMap(j,1),heatMap(j,2),heatMap(j,3),heatMap(j,4),heatMap(j,5)
    enddo
    close(1)
end program main




