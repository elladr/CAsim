PROGRAM CAFE

USE m_dengrowth
USE simParam

IMPLICIT NONE

CHARACTER(len=19) :: FI
!INTEGER:: iter


!iter=0
call compute_Param
call dendgrowth_init
call Output(iter)


!.................... main loop ...............

do iter=1, niter

call dendgrowth_main


 if (  ( mod(iter,iterOutFields) ==0 .and. (iter+1)<niter) .or. iter==niter )   call Output
 

end do

END PROGRAM CAFE

