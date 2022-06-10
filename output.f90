SUBROUTINE Output

use m_dengrowth
use simParam

implicit none


!integer, INTENT(INOUT):: iter
CHARACTER(len=19) :: FN
integer:: i,j, tDim, kk
integer:: px(xDim*yDim), py(xDim*yDim)
integer:: p_IvG(xDim*yDim), p_Iv(xDim*yDim)    


!	 ========================================================
!	 Write the components to a Tec. file 
!	 ========================================================

! ------------------------- Tecplot header for contours 




WRITE (FN, "(I8,A4)") iter, ".tec"
OPEN(250,FILE=FN)
write(250,1021)iter         
write(250,1022)xDim,yDim
 1021    format( 'TITLE = "Time = ',I8,'"',/, &
              'VARIABLES = "X", "Y" , "I_v_G", "I_v"') 
 1022    format( 'ZONE I=',i6,', J=',i6,', F=BLOCK')          
         
!   Tecplot data
    tDim=xDim*yDim
	
	

        do j=1, yDim
            do i=1, xDim
                kk=i + (j-1)*xDim
                px(kk)=i
                py(kk)=j

                p_IvG(kk)= I_v_G(i,j)
                p_Iv(kk)= I_v(i,j)
				
            end do
        end do
		
		
		


    
    write(250,310)(px(kk),kk=1,tDim)
    write(250,310)(py(kk),kk=1,tDim)  

    write(250,310)(p_IvG(kk),kk=1,tDim)
    write(250,310)(p_Iv(kk),kk=1,tDim)

    
  310    format(I6, ' ', I6, ' ', I6, ' ', I6, ' ', I6, ' ', I6, ' ', I6, ' ', I6)
! 310 format(1000(I6))


    close(250)



 
        
  RETURN      
 
      
END SUBROUTINE Output
