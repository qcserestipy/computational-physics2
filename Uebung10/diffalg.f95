program main
    use cp_mod

    implicit none
    integer, dimension(5,5) :: pointMap
    real, dimension(:), allocatable :: solutionVector
    real, dimension(:,:), allocatable :: matrix
    integer :: i,j
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

end program main




