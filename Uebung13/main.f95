program main 
    use cp_mod
    use FT_mod

    character(len=72) :: filename = "traj.xyz"
    type( system ), dimension(:), allocatable :: listofFrames
    real, dimension(:,:), allocatable :: corr

    call readtrjfile(filename, listofFrames,2220)
    call calcVelocities(listofFrames,2220)
    call autoCorrelationCalculator(listofFrames, 2220,corr)
end program main