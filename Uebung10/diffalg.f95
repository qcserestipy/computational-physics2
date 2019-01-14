program name

    interface
    subroutine write_matrix(a)
    real, dimension(:,:) :: a
    end subroutine write_matrix
    end interface


    
integer, parameter :: n = 4, m = 4

real, dimension(:), allocatable :: w_vec
real, dimension(1:3,1:3) :: pointmap
real, dimension(1:9,1:4) :: relationmap
integer, dimension (1:2) :: order2 = (/ 2, 1 /)

!setup matrix with indizes
!!use indextransformer
do i = 1, 3
    do j = 1, 3
        pointmap(i,j) = indexTransformer(i,j,n,m)
    enddo
enddo

call write_matrix(pointmap)

!setup relation matrix
do i = 1, 9
    do j = 1, 4
        relationmap(i,j) = 0.0
    enddo
enddo

do i = 1, 3
    do j = 1, 3
        if( (i+1) < 4 .and. (i+1) > 0) then
            relationmap(i,j) = pointmap(i+1,j)
        end if
    enddo
enddo



call write_matrix(relationmap)

!setup dlg for every point

!check neighbours if, value not on points put on right side (randbed)

!solve dlg

end program name

integer function indexTransformer(i,j,n,m)
integer :: i,j,n,m

indexTransformer = i + (m - 1 - j)*(n - 1)
return
end function indexTransformer


subroutine write_matrix(a)
    real, dimension(:,:) :: a
    write(*,*)
    
    do i = lbound(a,1), ubound(a,1)
       write(*,*) (a(i,j), j = lbound(a,2), ubound(a,2))
    end do
 end subroutine write_matrix
