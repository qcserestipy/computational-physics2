program main
    implicit none
    integer :: k,i,j
    
    call tuple2int(2,4,k)
    write(*,*) k

    call int2tuple(9,i,j)
    write(*,*) i,j
    
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

