module veloverlet_mod
    use cp_mod
    implicit none
    save

!##############################
!#    function definitions    #
!##############################
    contains


    real function calcDistance(atomA, atomB)
            implicit none
    
            type( atom ), intent(in) :: atomA, atomB
            real, dimension(3) :: distance
    
            distance(1) = distance(1) + atomA%coords%x - atomB%coords%x
            distance(2) = distance(2) + atomA%coords%y - atomB%coords%y
            distance(3) = distance(3) + atomA%coords%z - atomB%coords%z
            calcDistance = sqrt(distance(1)**2 + distance(2)**2 + distance(3)**2)
            return
    end function calcDistance

    subroutine timeloop(listofatoms, nmax, dt)
        implicit none
                type( atom ), dimension(:), allocatable, intent(inout) :: listofatoms
        real, INTENT(IN) :: dt
        real :: x
        integer, INTENT(IN) ::  nmax
        character(len=72) :: outfile = 'CO.trj' 
        character(len=72) :: energyFile = 'E.out'
        character(len=72) :: periodFile = 'period.out'
        character(len=72) :: velocityFile = 'velocity.out'
        real :: v = 0, f = 0.0d0, m = 15.9994, E = 0.0d0, t = 0
        real, dimension(:), allocatable :: distances
        integer :: nsteps, i, j

        call fileCleanUp(outfile)
        call fileCleanUp(energyFile)
        call fileCleanUp(periodFile)
        call fileCleanUp(velocityFile)

        write(*,*) "####################################"
        write(*,*) "#     Entering velocity verlet     #" 
        write(*,*) "####################################"
        nsteps = int(nmax/dt)
        allocate(distances(size(listofatoms)))
        x = calcDistance(listofatoms(1), listofatoms(2))
        call writexyzfile(outfile,listofatoms, size(listofatoms))
    
        do i = 1, size(distances)
            do j = i, size(distances)
                distances(i) = calcDistance(listofatoms(i),listofatoms(j))
            enddo
        enddo

        open(unit=3,file=energyFile)
        open(unit=4,file=periodFile)
        open(unit=5,file=velocityFile)
        do i = 1, nsteps
            v = v + f * dt / (2*m)
            x = x + v * dt
            f =  updateForce(x)
            E = updateEnergy(x)
            t = t + dt
            listofatoms(2)%coords%z = x
            call writetrjfile(outfile,listofatoms, size(listofatoms))
            write(3,*) x, E
            write(4,*) t, x
            write(5,*) t, v
        enddo
        close(5)
        close(4)
        close(3)
    end subroutine timeloop
    
    real function updateForce(x)
        implicit none
    real :: x
    real, parameter :: D = 0.493172, a = 2.05666, R_0 = 1.20106
        updateForce = - 2 * a * D * exp(a*(R_0 - x)) * (1 - exp(a*(R_0 - x)))
    end function updateForce
    
    
    real function updateEnergy(x)
        implicit none
    real :: x
    real, parameter :: D = 0.493172, a = 2.05666, R_0 = 1.20106
        updateEnergy = D*(( 1 - exp(-a*(x - R_0)) )**2)
    end function updateEnergy

end module veloverlet_mod
    