module FT_mod
    use cp_mod
    implicit none
    save


    !##############################
    !#    function definitions    #
    !##############################
    contains

    subroutine discreteFourierTransformer(f_x, f_k, N)
        implicit none
        integer, intent(in) :: N
        complex*16, dimension(N), intent(in) :: f_x
        complex*16, dimension(N), intent(out) :: f_k
        real*16, parameter :: PI = 4 * ATAN(1.)
        integer :: i,j

        write(*,*) "####################################"
        write(*,*) "#       Entering F-transform       #" 
        write(*,*) "####################################"

        do i=1,N
            f_k(i) = 0
            do j=1,N
                f_k(i) = f_k(i) + f_x(j) * cexp(cmplx(0, -2*PI*i*j/N))
            end do
            !write(*,*) i, f_k(i)
        end do
    end subroutine discreteFourierTransformer




    subroutine powerSpectrum(dataPoints,spec_points,outfile)
        implicit none
        character(len=72), intent(in) :: outfile
        integer, intent(in) :: spec_points
        real, dimension(spec_points), intent(in) :: dataPoints 
        complex*16, dimension(:), allocatable :: f_x
        complex*16, dimension(:), allocatable :: f_k
        real*16, parameter :: PI = 4 * ATAN(1.)
        real, dimension(:), allocatable :: power_spec, omega
        real :: dt = 0.002
        integer :: i

        call fileCleanUp(outfile)

        allocate(omega(spec_points))
        allocate(power_spec(spec_points))
        allocate(f_k(spec_points))
        allocate(f_x(spec_points))
        

        do i = 1, spec_points
            f_x(i) = cmplx(0,dataPoints(i))
        enddo

        call discreteFourierTransformer(f_x, f_k, spec_points)

        do i = 1, spec_points
            omega(i) = i * 0.98652!i * ((2*PI)/(spec_points*dt))
        enddo

        write(*,*) "#########################################"
        write(*,*) "#        Entering power spectrum        #" 
        write(*,*) "#########################################"
        open(unit=1,file=outfile)
        do i = 1, spec_points
            power_spec(i) = abs(f_k(i))**2
            write(1,*) omega(i), power_spec(i)
        enddo
        close(1)
    end subroutine powerSpectrum




    subroutine autoCorrelationCalculator(listofFrames, nframes,correlation)
        implicit none
        type( system ), dimension(:), allocatable, intent(in) :: listofFrames
        real, dimension(:,:), allocatable, intent(out) :: correlation
        integer, intent(in) :: nframes
        integer :: i,j,k, natoms

        natoms = size(listofFrames(1)%subsystem)
        allocate(correlation(nframes,natoms))

        do i = 1, nframes - 1
            do j = 1, natoms
                k = i + 1
                correlation(i,j) = (real(listofFrames(i)%subsystem(j)%s_velo) &
                                  * real(listofFrames(k)%subsystem(j)%s_velo) )
            enddo
        enddo
    end subroutine autoCorrelationCalculator




    subroutine greenKuboCalculator(listofFrames, nframes,diff)
        implicit none
        type( system ), dimension(:), allocatable, intent(in) :: listofFrames
        real, intent(out) :: diff
        integer, intent(in) :: nframes
        integer :: i,j,k,natoms

        natoms = size(listofFrames(1)%subsystem)
        

        do i = 1, nframes - 1
            do j = 1, natoms
                k = i + 1
                diff = diff + (real(listofFrames(k)%subsystem(j)%s_velo) &
                            *  real(listofFrames(i)%subsystem(j)%s_velo) )
            enddo
        enddo
        write(*,*) "Diffusionskoeffizient Green-Kubo: ", diff / (3.0*natoms*nframes)
    end subroutine greenKuboCalculator




    subroutine einsteinCalculator(listofFrames, nframes)
        implicit none
        type( system ), dimension(:), allocatable, intent(in) :: listofFrames
        integer, intent(in) :: nframes
        real :: sum = 0.0d0
        integer :: i,j,k, natoms

        natoms = size(listofFrames(1)%subsystem)

        open(unit=1,file="einstein.out")
        do i = 1, nframes -1 
            do j = 1, natoms
                k = i + 1
                sum = sum + distance(listofFrames(i)%subsystem(j),listofFrames(k)%subsystem(j))**2
            enddo
            write(1,*) i,",", sum / real(nframes)
        enddo
        close(1)    
    end subroutine einsteinCalculator




    subroutine radialDistribution(listofFrames,nframes,boxL,corr)
        implicit none
        type( system ), dimension(:), allocatable, intent(in) :: listofFrames
        integer, intent(in) :: nframes
        real, INTENT(IN) :: boxL
        
        real :: dR = 0.49326
        real*4, DIMENSION(:,:,:), allocatable :: hitlist
        type( system ), dimension(:), allocatable :: listofOxygens
        type( atom ), dimension(:), allocatable :: listofatoms
        integer :: i,j,k,natoms,l
        real :: disti, disto, distAB, vi, vo
        real :: sum = 0.0d0
        real, DIMENSION(:,:), ALLOCATABLE :: local 
        real, DIMENSION(:), ALLOCATABLE :: time_local 
        real, DIMENSION(:), ALLOCATABLE,intent(out) :: corr


        natoms = size(listofFrames(1)%subsystem)
        allocate(hitlist(2220,32,20))
        allocate(listofOxygens(nframes))
        allocate(listofatoms(32))
        hitlist = 0

        !collect all O atoms in one system list
        do i = 1, nframes
            do j = 1, natoms
                if (listofFrames(i)%subsystem(j)%name == 'O') then
                    listofatoms(j) = listofFrames(i)%subsystem(j)
                end if
            enddo
            listofOxygens(i)%subsystem = listofatoms
        enddo

        !count all O atoms for every frame within the discrete box volumes
        do l = 1, 2220     
            do i = 1, 32
                do k = 1, 20
                    do j = 1, 32
                        if (i == j) then
                            continue
                        else
                            disti = (k-1)*dR
                            disto = k*dR
                            distAB = periodicDistance(listofOxygens(l)%subsystem(i), listofOxygens(l)%subsystem(j), boxL)
                            if (distAB < disto .and. distAB > disti) then
                                hitlist(l,i,k) = hitlist(l,i,k) + 1
                            end if
                        end if
                    enddo
                enddo
            enddo
        enddo
        
        !scale every count with discrete volume, still contains issues
        !do l = 1, 2220
        !    do i = 1, 32
        !        do k = 1, 10
        !            disti = (k-1)*dR
        !            disto = k*dR
        !            vi = vol(disti)
        !            vo = vol(disto)
        !            hitlist(l,i,k) = hitlist(l,i,k)
        !        enddo
        !    enddo
        !enddo

        !average over all particles and average over all time steps
        allocate(local(nframes,20))
        allocate(time_local(20))

        do l = 1, 2220
            do k = 1, 20
                disti = (k-1)*dR
                    disto = k*dR
                    vi = vol(disti)
                    vo = vol(disto)
                do i = 1, 32
                    sum = sum + hitlist(l,i,k)
                enddo
                sum = sum / 32.0
                local(l,k) = sum
                sum = 0.0d0
            enddo
        enddo


        do k = 1,20
            do l = 1, 2220
                 sum = sum + local(l,k)
            enddo
            sum = sum / 2220.0
            time_local(k) = sum
            sum = 0.0d0
        enddo

        !calculate correlation of density
        allocate(corr(19))
        do k = 1, 19
            j = k+1
            corr(k) = time_local(k) * time_local(j)
        enddo

    
        do k = 1,20
            sum = sum + time_local(k)
        enddo

        !io
        open(unit=1,file="corr.out")
        do k = 1,19
            corr(k) = corr(k) / sum
            write(1,*) k*dR, corr(k)
        enddo
        close(1)

    end subroutine radialDistribution


    real function vol(r)
        implicit none
        real :: r
        real*16, parameter :: PI = 4 * ATAN(1.)

        vol = (4.0/3.0)*PI*(r**3)
    end function vol



!    subroutine gr(switch)
!        implicit none
!
!
!        if (switch .eq. 0) then
!            ngr = 0
!            delg=box/(2*nhis)
!            do i=0,nhis
!                g(i) = 0
!            enddo
!        else if (switch .eq. 1) then
!            ngr=ngr+1
!            do i = 1, npart - 1
!                do j =i+1,npart
!                    xr=x(i)-x(j)
!                    xr=xr-box*nint(xr/box)
!                    r=sqrt(xr**2)
!                    if (r .lt. box/2) then
!                        ig = int(r/delg)
!                        g(ig)=g(ig)+2
!                    endif
!                enddo
!            enddo
!        else if (switch .eq. 2) then
!            do i = 1, nhis
!                r = delg*(i+0.5)
!                vb=((i+1)**3-i**3)*delg**3
!                nid=(4/3)*PI*vb*rho
!                g(i)=g(i)/(ngr*npart*nid)
!            enddo
!        endif
!        RETURN
!        end subroutine gr
end module FT_mod