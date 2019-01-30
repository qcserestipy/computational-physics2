program main 
    use cp_mod
    use FT_mod

    character(len=72) :: filename = "traj.xyz"
    character(len=72) :: outfile = "powerSpec.out"
    type( system ), dimension(:), allocatable :: listofFrames
    real, dimension(:,:), allocatable :: corr

    call readtrjfile(filename, listofFrames,2220)
    call calcVelocities(listofFrames,2220)
    call autoCorrelationCalculator(listofFrames, 2220,corr)
    call powerSpectrum(corr,2220,outfile)
end program main