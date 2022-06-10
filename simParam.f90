MODULE simParam

implicit none
INTEGER:: xDim, yDim, Y_nucleation

REAL:: D, Gamma, kcoeff, ml, c0
REAL:: Acoef, ncoef
REAL:: GradT, Vpull 
REAL:: TInit, UNDERCOOL_0, TLequi, TLpure
REAL:: FracGrain, ANGLE1DEG, ANGLE2DEG  
INTEGER:: iGB1, iGB2, nlayer
REAL:: TotalTime, Kdt, dt, dx
INTEGER:: NwriteF, niter, iterOutFields, iter
REAL, DIMENSION(2) :: Alpha
REAL, PARAMETER :: PI = 3.1415927
LOGICAL, PARAMETER:: Temp_CC= .TRUE.        !if true use temp at cell center else use temp at Growth center
LOGICAL, PARAMETER:: Single_Crys= .TRUE.        !if true model a singel equiaxed grain
!integer:: iter



CONTAINS

SUBROUTINE compute_Param




	
!.......................................Physical properties
c0= 1.3      ! UnitCompo
ml= -3.02     ! K/UnitCompo
kcoeff= 0.1
Gamma= 6.4e-8	! GIBBSTHOMSON K.m
!Eps4= 0.014   ! ANISOTROPY	
D= 1270.   ! micron^2/s
ncoef= 2
Acoef= 100        !1.08*1.e-4*1.e6   ! micron / (s K^n)

! ........................................Process
GradT= 1.*250.*1.e-6    ! K/micron
Vpull= 400.         ! micron/s
UNDERCOOL_0= 2.
TLpure= 331.24                    ! Tlpure is the melting temperature of pure metal
                          ! Tunder is under cooling
TLequi = TLpure + ml*c0           ! TLequi is the temperature at equalibrium at solute 
TInit = TLequi - UNDERCOOL_0           ! TInit is the initial temperature 

!.......................................Computational condition
iter=1
TotalTime=40.      ! s
Kdt=0.1
dx=100.    !micron
dt= 0.025     !Kdt*dx/Vpull	
NwriteF= nint(TotalTime)	
niter=int(TotalTime/dt)
iterOutFields=niter/NwriteF	
xDim= 100
yDim= 100
Y_nucleation= 50


!........................................Grain
FracGrain= 1.
ANGLE1DEG= 30   !90-30.
ANGLE2DEG= 30  !90-10.
Alpha=	(/ ANGLE1DEG*PI/180.,ANGLE2DEG*PI/180. /)	! rotation angle of anisotropy
iGB1= nint((xDim+1)*FracGrain*.5)-1
iGB2= nint((xDim+1)*(1.-FracGrain*.5))-1
nlayer=1




open(1, file="param.txt")
write(1,*) "niter=", niter
write(1,*) "iGB1=", iGB1
write(1,*) "iGB2=", iGB2
write(1,*) "nlayer=", nlayer
!close(1)



END SUBROUTINE compute_Param

END MODULE simParam






