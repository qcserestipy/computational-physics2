program main 
    use cp_mod
    use FT_mod

    character(len=72) :: filename = 'velos.txt'
    character(len=72) :: outA ="power_spec_A.out"
    character(len=72) :: outB ="power_spec_B.out"
    character(len=72) :: outC ="power_spec_C.out"
    character(len=72) :: outD ="power_spec_D.out"

    real, dimension(18000) :: t_velo
    real, dimension(2000) :: corrA
    real, dimension(2000) :: corrB
    real, dimension(2000) :: corrC
    real, dimension(2000) :: corrG

    call MDReader(filename,t_velo)
    call autoCorrelationCalculator(t_velo,corrA,corrB,corrC,corrG)

    call fileCleanUp(outA)
    call fileCleanUp(outB)
    call fileCleanUp(outC)
    call fileCleanUp(outD)
    call powerSpectrum(corrA,2000,outA)
    call powerSpectrum(corrB,2000,outB)
    call powerSpectrum(corrC,2000,outC)
    call powerSpectrum(corrG,2000,outD)

end program main