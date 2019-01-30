program main 
    use cp_mod
    use FT_mod

    character(len=72) :: filename = "traj.xyz"
    character(len=72) :: outfile = "powerSpec.out"
    type( system ), dimension(:), allocatable :: listofFrames
    real, dimension(:,:), allocatable :: corr
    real :: diffusion = 0.0d0

    call readtrjfile(filename, listofFrames,2220)
    call calcVelocities(listofFrames,2220)
    call autoCorrelationCalculator(listofFrames, 2220,corr)
    call powerSpectrum(corr,2220,outfile)
    call greenKuboCalculator(listofFrames, 2220,diffusion)

    call einsteinCalculator(listofFrames, 2220)
end program main