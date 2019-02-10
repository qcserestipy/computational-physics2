program main 
    use cp_mod
    use FT_mod

    character(len=72) :: filename = "traj.xyz"
    character(len=72) :: outfile = "powerSpec.out"
    type( system ), dimension(:), allocatable :: listofFrames
    real, dimension(:), allocatable :: corr
    real :: diffusion = 0.0d0
    real :: boxL = 9.8652

    call readtrjfile(filename, listofFrames,2220)
    !call calcVelocities(listofFrames,2220)
    !call autoCorrelationCalculator(listofFrames, 2220,corr)
    
    !call greenKuboCalculator(listofFrames, 2220,diffusion)

    !call einsteinCalculator(listofFrames, 2220)
    call radialDistribution(listofFrames,2220,boxL,corr)
    !call powerSpectrum(corr,10,outfile)

end program main