program test
use cp_mod

implicit none

    character(len=72) :: infile = 'water.xyz'
    character(len=72) :: outfile = 'output.xyz'
    type ( atom ), dimension(:), allocatable :: water

    call readxyzfile(infile, water)
    call writexyzfile(outfile, water, size(water))

end program test