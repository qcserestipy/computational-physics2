program test

real, dimension(2,2) :: arr(2,2)

arr = reshape ((/2,2,1,2/),(/2,2/))

write(*,*) arr()
end program test
