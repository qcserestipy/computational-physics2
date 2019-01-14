!================================================
PROGRAM mc
!================================================

INTEGER :: n, nacc = 0, nstep = 300000
REAL :: e = 0 
REAL, dimension(:), allocatable :: x
REAL, dimension(:), allocatable :: y
REAL, dimension(:), allocatable :: z
CHARACTER(LEN=2), dimension(:), allocatable :: name

!read system 

open(unit=10,file='au13-ini.xyz', iostat=io_error)

      if (io_error == 0) then
        read(10,*) n
        allocate (x(n))
        allocate (y(n))
        allocate (z(n))
        allocate (name(n))

        do k = 1, n
             read(10,*) name(k),  x(k), y(k), z(k)
        end do

        write(*,*) "Reading file, getting atoms with following coordinates:"

        write(*,*) n
        do k = 1, n
             write(*,*) name(k),  x(k), y(k), z(k)
        end do

        write(*,*) "Reading completed!"

       else
           write(*,*) 'Beim OEffenen der Datei ist ein Fehler Nr.', &
                       io_error,' aufgetreten'
       end if
close(unit=10)





write(*,*) "nstep", "nacc", "gamma" !gamma = Akzeptanzrate
write(*,*) "=============================================="


! nstep MC steps 
DO istep=1,nstep 


! randomly displace atoms
   CALL mcmove(x,y,z,n,nacc,e)

   write(*,*) "istep" ,istep, "nacc", nacc, "ratio", REAL(nacc) / REAL(istep)   








! write output
open (unit=20,file="trajectory.xyz",position='append')
  write(20,*) n
  write(20,*)
  do k = 1, n
    write(20,*) name(k), x(k), y(k), z(k)
  end do
close(20)

ENDDO

write(*,*) "Program ended. All clean and tidy!"

!================================================
END PROGRAM mc
!================================================


!================================================
SUBROUTINE mcmove(x,y,z,n,nacc,e)
!================================================
INTEGER :: n,nacc 
REAL :: e, eold , delta = 0.13, t = 2000, boltz = 0
REAL , PARAMETER :: kb = 1.38d-23
REAL, dimension(n) :: x
REAL, dimension(n) :: y
REAL, dimension(n) :: z
REAL, dimension(n) :: xnew
REAL, dimension(n) :: ynew
REAL, dimension(n) :: znew



! calculate energy at old positions
CALL energy(x,y,z,n,e)
eold=e

! random displacement of atoms: 
DO i=1,n

   CALL random_number(rand)
   xnew(i)=x(i)+delta*(rand-0.5d0)
   CALL random_number(rand)
   ynew(i)=y(i)+delta*(rand-0.5d0)
   CALL random_number(rand)
   znew(i)=z(i)+delta*(rand-0.5d0)

ENDDO

! calculate energy at new positions
CALL energy(xnew,ynew,znew,n,e)



! calculate Boltzmann factor: 
! kb = Boltzmann constant, t = temperature

boltz=EXP(-(e-eold)/(kb*t))
!write(*,*) "boltz", boltz, "rand", rand, "e" ,e, "eold", eold


! accept or reject new positions
CALL random_number(rand)
IF (rand.LT.boltz) THEN

   DO i=1,n

      x(i)=xnew(i)
      y(i)=ynew(i)
      z(i)=znew(i)

   ENDDO

! count accepted trial moves
   nacc=nacc+1

ELSE

   e=eold

ENDIF

RETURN
!================================================
END SUBROUTINE mcmove
!================================================

!================================================
SUBROUTINE energy(x,y,z,n,e)
!================================================

!IMPLICIT NONE
INTEGER :: n,i,j
REAL :: e,e2 = 0,r,dx,dy,dz
REAL, PARAMETER :: a=4.08,c=34.408,eps=2.04968D-21 !IN JOULE
!REAL, dimension(n) :: x
!REAL, dimension(n) :: y
!REAL, dimension(n) :: z
REAL :: x(n), y(n), z(n)

e=0.0d0

DO i=1,n

   e2=0.0d0

   DO j=1,n

      IF(i.NE.j)THEN
         dx=x(i)-x(j)
         dy=y(i)-y(j)
         dz=z(i)-z(j)

         r=SQRT(dx**2+dy**2+dz**2)

         e=e+0.5d0*(a/r)**10

         e2=e2+(a/r)**8

      ENDIF
   ENDDO
   e=e-c*SQRT(e2)
ENDDO
! Umrechnen der Energie in Joule
e=e*eps
RETURN
!================================================
END SUBROUTINE energy
!================================================
