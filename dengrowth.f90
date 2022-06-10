!---------------------------------------------------------------------
SUBROUTINE dendgrowth_init
!---------------------------------------------------------------------
  
    USE m_dengrowth
    USE SimParam

  
    IMPLICIT NONE
    INTEGER :: i, j, ii, jj
    INTEGER:: nuclei_y, nuclei_x1G1, nuclei_x2G1, nuclei_xG2

    ALLOCATE(I_v(xDim, yDim))
    ALLOCATE(I_v_G(xDim,yDim))
    ALLOCATE(width_cell_old(xDim,yDim))
    ALLOCATE(width_cell(xDim,yDim))
    ALLOCATE(V_tip(xDim,yDim))
    ALLOCATE(tp_G(xDim,yDim))
    ALLOCATE(tp(yDim))
    ALLOCATE(tp_init(yDim))
    ALLOCATE(Teta(xDim,yDim))
    ALLOCATE(x_G(xDim,yDim))
    ALLOCATE(y_G(xDim,yDim))
    ALLOCATE(x_S10(xDim,yDim))
    ALLOCATE(y_S10(xDim,yDim))
    ALLOCATE(x_S01(xDim,yDim))
    ALLOCATE(y_S01(xDim,yDim))
    ALLOCATE(x_Sn10(xDim,yDim))
    ALLOCATE(y_Sn10(xDim,yDim))
    ALLOCATE(x_S0n1(xDim,yDim))
    ALLOCATE(y_S0n1(xDim,yDim))


    I_v(:,:)= 0     !liquid=0, solid=1
    I_v_G(:,:)=0

    width_cell_old(:,:)= 0.
    width_cell(:,:)= 0.

  !  V_tip_old(:,:)= 0.
    V_tip(:,:)= 0.



    tp_init(:)= 0.
    tp(:)= 0.
    tp_G(:,:)= 0.


    DO j = 1, yDim
        
        tp_init(j) = TInit + GradT * (j - Y_nucleation) * dx  
        tp(j) = tp_init(j)                         !temp at the cell center

        tp_G(:,j) = TInit + GradT * (j - Y_nucleation) * dx     !temp at the Growth Center
     
    END DO




    Teta(:,:)= 0.

    x_G(:,:)= 0.
    y_G(:,:)= 0.
    x_S10(:,:)= 0.
    y_S10(:,:)= 0.
    x_S01(:,:)= 0.
    y_S01(:,:)= 0.
    x_Sn10(:,:)= 0.
    y_Sn10(:,:)= 0.
    x_S0n1(:,:)= 0.
    y_S0n1(:,:)= 0.






    nuclei_y = Y_nucleation  


    IF (Single_Crys)     THEN


        nuclei_x1G1 = iGB1
        nuclei_x2G1= iGB1
        nuclei_xG2= iGB1 

        ELSE

        nuclei_x1G1 = iGB1/2
        nuclei_x2G1= iGB2+iGB1/2 
        nuclei_xG2= (iGB1+iGB2)/2 

    ENDIF 

    write(1,*) "nuclei_y=", nuclei_y
    write(1,*) "nuclei_x1G1=", nuclei_x1G1
    write(1,*) "nuclei_x2G1=", nuclei_x2G1
    write(1,*) "nuclei_xG2", nuclei_xG2
    close(1)            

    I_v(nuclei_x1G1, nuclei_y) = 1
    I_v(nuclei_x2G1, nuclei_y) = 1
    I_v(nuclei_xG2, nuclei_y) = 1

    x_G(nuclei_x1G1, nuclei_y) = nuclei_x1G1 * dx
    y_G(nuclei_x1G1, nuclei_y) = nuclei_y * dx

    x_G(nuclei_x2G1,nuclei_y) = nuclei_x2G1 * dx
    y_G(nuclei_x2G1,nuclei_y) = nuclei_y * dx   
    
    x_G(nuclei_xG2,nuclei_y) = nuclei_xG2 * dx
    y_G(nuclei_xG2,nuclei_y) = nuclei_y * dx





   

    DO j=1,yDim
        DO i=1,xDim
            IF ( I_v(i,j) == 1 ) THEN
        
                IF (i<iGB1 .or. i>=iGB2) THEN
                    I_v_G(i,j)=1
                    Teta(i,j)= Alpha(1)


                ELSE
                    I_v_G(i,j)=2
                    Teta(i,j)= Alpha(2)
                ENDIF 

                print *, "i & j solids=", i,j
                print *, "tp_init j=", tp_init(j)
                print *, "TInit=", TLequi
                print *, "2nd term=", ( TLequi- tp_init(j) )
                print *, "2nd term2=", ( TLequi- tp_G(i,j) )  


                IF (Temp_CC ) THEN

                    V_tip(i,j) = Acoef* ( TLequi- tp(j) )** ncoef

                    ELSE
                
                    V_tip(i,j) = Acoef* ( TLequi- tp_G(i,j) )** ncoef

                ENDIF

                width_cell_old(i,j)= ( 2./sqrt(2.) )*  V_tip(i,j) *dt 
                width_cell(i,j)= width_cell_old(i,j)

                x_S10(i,j)= x_G(i,j)+0.5*sqrt(2.)*width_cell(i,j)*cos(Teta(i,j))
                y_S10(i,j)= y_G(i,j)+0.5*sqrt(2.)*width_cell(i,j)*sin(Teta(i,j))
                x_S01(i,j)= x_G(i,j)-0.5*sqrt(2.)*width_cell(i,j)*sin(Teta(i,j))


                print *, "i & j=", i, j
                print *, "V_tip_nuclei=", V_tip(i,j)
                print *, "width_cell_nuclei=", width_cell(i,j)
                print *, "xG=", x_G(i,j)
                print *, "yG=", y_G(i,j)
                print *, "xS10=", x_S10(i,j)



            ENDIF
        
        END DO
    END DO


 
END SUBROUTINE dendgrowth_init



!=========================================================================================================
!!----------calculate the dendrite growth
!=========================================================================================================
SUBROUTINE dendgrowth_main

 USE m_dengrowth
 USE SimParam, ONLY: xDim, yDim, dt, iter, Vpull, GradT


 IMPLICIT NONE
 INTEGER :: i, j






CALL dendgrowth_CA

DO j = 1, yDim               
    tp(j) =  tp_init(j)  - Vpull*GradT *(iter-0.)* dt          
END DO

CALL cell_growth_kinetic

width_cell_old = width_cell 

    
END SUBROUTINE dendgrowth_main




!=========================================================================================================
!!----------calculate cell width
!=========================================================================================================
SUBROUTINE cell_growth_kinetic

 USE m_dengrowth
 USE SimParam


 IMPLICIT NONE
 INTEGER :: i, j, ns, ii, jj

    ns=0


    DO j=1,yDim
        DO i=1,xDim
            IF ( I_v(i,j) == 1 ) THEN

                ns = 0  

                DO jj = -nlayer, nlayer
                    DO ii = -nlayer, nlayer
                        IF(I_v(i+ii,j+jj)==1  ) THEN
                            ns= ns+1
                        ENDIF
                    ENDDO
                ENDDO

                IF (ns .lt. 9) THEN


                IF (Temp_CC ) THEN

                    V_tip(i,j) = Acoef* ( TLequi- tp(j) )** ncoef

                    ELSE

                    tp_G(i,j) =  TInit + GradT * (y_G(i,j)- 2.*dx) -  Vpull*GradT *(iter-0)* dt
                    V_tip(i,j) = Acoef* ( TLequi- tp_G(i,j) )** ncoef

                ENDIF

                width_cell(i,j)= width_cell_old(i,j) + ( 2./sqrt(2.) )*  V_tip(i,j) *dt 

                x_S10(i,j)= x_G(i,j)+0.5*sqrt(2.)*width_cell(i,j)*cos(Teta(i,j))
                y_S10(i,j)= y_G(i,j)+0.5*sqrt(2.)*width_cell(i,j)*sin(Teta(i,j))
                x_S01(i,j)= x_G(i,j)-0.5*sqrt(2.)*width_cell(i,j)*sin(Teta(i,j))


!                print *, "i & j=", i, j
!                print *, "V_tip_nuclei=", V_tip(i,j)
!                print *, "width_cell_nuclei=", width_cell(i,j)
!                print *, "xG=", x_G(i,j)
!                print *, "yG=", y_G(i,j)
!                print *, "xS10=", x_S10(i,j)

                ENDIF



            ENDIF
        
        END DO
    END DO



END SUBROUTINE cell_growth_kinetic



!    ========================================================
!    Calculation the dendrite growth
!    ========================================================
SUBROUTINE dendgrowth_CA
USE m_dengrowth
USE simParam
IMPLICIT NONE

REAL, DIMENSION(5) :: coord_S_help
REAL, DIMENSION(10):: G_and_S_coords
INTEGER, DIMENSION(8):: solid_ngbr_x, solid_ngbr_y
REAL, DIMENSION(8):: width_compare
INTEGER:: ii, jj, i, j, ns, k, index
REAL:: xS_nearest, yS_nearest, nearest_vertice_no
REAL:: width_cell_0
REAL:: Teta_help, Teta_C_M
REAL:: w_G_C, w_Cprime_G
REAL:: L1_m, L2_m


coord_S_help= 0.
G_and_S_coords = 0.
nearest_vertice_no= 0.
xS_nearest= 0.
yS_nearest= 0.
width_cell_0= 0.
Teta_help= 0.
Teta_C_M= 0.
w_Cprime_G= 0.
w_G_C= 0.
L1_m= 0.
L2_m= 0.


solid_ngbr_x(:)= 0
solid_ngbr_y(:)= 0
width_compare(:)= 0.

ns=0 


      DO j = 1+nlayer, yDim-nlayer
        DO i = 1+nlayer, xDim-nlayer

            IF(I_v(i,j)==0) THEN
                ns = 0  
                solid_ngbr_x(:)= 0
                solid_ngbr_y(:)= 0
                width_compare(:)= 0. 

                DO jj = -nlayer, nlayer
                    DO ii = -nlayer, nlayer
                        IF(I_v(i+ii,j+jj)==1  .and. I_v(i,j)==0) THEN

                        coord_S_help= nearest_vertice_to_liquid_cell_cntr(i,j,i+ii,j+jj)
                        xS_nearest= coord_S_help(1)
                        yS_nearest= coord_S_help(2)
                        nearest_vertice_no= coord_S_help(3)



                        Teta_help=angle_between_liquid_cell_cntr_and_nearest_vertice(i,j, i+ii, j+jj, &
                                                                                xS_nearest,yS_nearest)

                        Teta_C_M= 45*PI/180- Teta_help

                     !   w_Cprime_G= 0.5* width_cell(i+ii, j+jj) / cos(Teta_C_M)
                        w_Cprime_G= 0.5* width_cell(i+ii, j+jj) 
                        w_G_C= cos(Teta_C_M)* sqrt ( (i*dx -x_G(i+ii,j+jj)) **2. + (j*dx-y_G(i+ii,j+jj)) **2. )

    !                    print *, "Teta_C_M=", Teta_C_M
    !                    print *, "w_Cprime_G=", w_Cprime_G
    !                    print *, "w_G_C", w_G_C

                        IF ( w_Cprime_G .ge. w_G_C ) THEN

                            ns = ns + 1
                            solid_ngbr_x(ns)=ii
                            solid_ngbr_y(ns)=jj
                            width_compare(ns)= w_Cprime_G - w_G_C
                        ENDIF

                        ENDIF

                    ENDDO
                ENDDO




                IF (ns .ge. 1) THEN

                    DO k=1,ns
                        ii= solid_ngbr_x(k)
                        jj= solid_ngbr_y(k)
                        
                    ENDDO



                    ! print "maxloc", maxloc(width_compare, DIM=1)
                    index= maxloc(width_compare, DIM=1)
                    ii= solid_ngbr_x(index)
                    jj= solid_ngbr_y(index)
                    print *, "widthcompare=", width_compare(index) 

                  
                    I_v(i,j)= I_v(i+ii,j+jj)
                    I_v_G(i,j)= I_v_G(i+ii,j+jj) 
                    Teta(i,j)= Teta(i+ii,j+jj)
                    IF ( I_v(i,j) .eq. 0 ) THEN
                        print *, "I_v wrong", I_v(i,j) 
                        pause
                    ENDIF

!                   IF (ns .ge. 2) THEN
!                       print *, "iter=", iter
!                       print *, "has", ns,  "solid neghbr"
!                       print *, "solid ngbr x=", solid_ngbr_x
!                       print *, "solid ngbr y=", solid_ngbr_y
!
!                   ENDIF


                    coord_S_help= nearest_vertice_to_liquid_cell_cntr(i,j,i+ii,j+jj)
                    xS_nearest= coord_S_help(1)
                    yS_nearest= coord_S_help(2)
                    nearest_vertice_no= coord_S_help(3)


                    L1_m= coord_S_help(4)
                    L2_m= coord_S_help(5)
                    width_cell_0= calculate_cell_width(L1_m, L2_m)
                    width_cell(i,j)= width_cell_0
                    width_cell_old(i,j)= width_cell_0
                    print *, "widthcell0=", width_cell_0 
                
                    G_and_S_coords= calculate_cntr_n_vertices(nearest_vertice_no, xS_nearest, yS_nearest, &
                                                                width_cell_0, Teta(i,j))
                    
                    x_G(i,j)= G_and_S_coords(1)
                    y_G(i,j)= G_and_S_coords(2)
                    x_S10(i,j)= G_and_S_coords(3)
                    y_S10(i,j)= G_and_S_coords(4)


                ENDIF


            
            ENDIF
       
        END DO 
    END DO 



CONTAINS



REAL FUNCTION calculate_cell_width(L1_min, L2_min)

USE simParam
USE m_dengrowth
IMPLICIT NONE 
REAL, INTENT(IN):: L1_min, L2_min



        calculate_cell_width= min ( L1_min, dx*sqrt(2.) ) + min ( L2_min, dx*sqrt(2.) )



END FUNCTION calculate_cell_width





FUNCTION calculate_cntr_n_vertices(n_vert_no, xS_near, yS_near, w_cell, teta_n)
                   
    IMPLICIT NONE
 !   INTEGER, INTENT(IN):: x, y
    REAL, INTENT(IN):: n_vert_no, xS_near, yS_near, w_cell, teta_n
    REAL, DIMENSION(10):: calculate_cntr_n_vertices
  !  INTEGER:: i,j
    REAL:: xx_s10, yy_s10, xx_s01, yy_s01, xx_sn10, yy_sn10, xx_s0n1, yy_s0n1, xx_G, yy_G


    IF (n_vert_no==1.)  THEN
        xx_s10= xS_near
        yy_s10= yS_near
        xx_sn10= xS_near- sqrt(2.) * w_cell* cos(teta_n)
        yy_sn10= yS_near- sqrt(2.) * w_cell* sin(teta_n)
        xx_G= xx_s10- 0.5*sqrt(2.) * w_cell* cos(teta_n)


    ELSEIF (n_vert_no==2.) THEN
        xx_s01= xS_near
        yy_s01= yS_near
        xx_s0n1= xS_near+ sqrt(2.) * w_cell* sin(teta_n)
        yy_s0n1= yS_near- sqrt(2.) * w_cell* cos(teta_n)



     ELSEIF (n_vert_no==3.)  THEN
        xx_sn10= xS_near
        yy_sn10= yS_near
        xx_s10= xS_near+ sqrt(2.) * w_cell* cos(teta_n)
        yy_s10= yS_near+ sqrt(2.) * w_cell* sin(teta_n)
        xx_G= xx_sn10+ 0.5*sqrt(2.) * w_cell* cos(teta_n)



    ELSEIF (n_vert_no==4.) THEN
        xx_s0n1= xS_near
        yy_s0n1= yS_near
        xx_s01= xS_near- sqrt(2.) * w_cell* sin(teta_n)
        yy_s01= yS_near+ sqrt(2.) * w_cell* cos(teta_n)
        xx_G= xx_s0n1- 0.5*sqrt(2.) * w_cell* sin(teta_n)


    ENDIF


    calculate_cntr_n_vertices(1)= xx_G
    calculate_cntr_n_vertices(2)= yy_G
    calculate_cntr_n_vertices(3)= xx_s10
    calculate_cntr_n_vertices(4)= yy_s10
    calculate_cntr_n_vertices(5)= xx_s01
    calculate_cntr_n_vertices(6)= yy_s01
    calculate_cntr_n_vertices(7)= xx_sn10



END FUNCTION calculate_cntr_n_vertices





FUNCTION nearest_vertice_to_liquid_cell_cntr(x,y,xx,yy)


    !USE m_dengrowth, ONLY :  x_S10, y_S10, x_S01, y_S01, x_Sn10, y_Sn10, x_S0n1, y_S0n1

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: y, x, xx, yy   
    REAL:: d1, d2, d3, d4
    REAL, DIMENSION(5) :: nearest_vertice_to_liquid_cell_cntr
   ! INTEGER :: i,j


    d1 = sqrt(     (x*dx-x_S10(xx,yy)) **2.   +   (y*dx-y_S10(xx,yy)) **2.     )
    d2 = sqrt(     (x*dx-x_S01(xx,yy)) **2.   +   (y*dx-y_S01(xx,yy)) **2.     )
    d3 = sqrt(     (x*dx-x_Sn10(xx,yy)) **2.   +   (y*dx-y_Sn10(xx,yy)) **2.     )
    d4 = sqrt(     (x*dx-x_S0n1(xx,yy)) **2.   +   (y*dx-y_S0n1(xx,yy)) **2.     )

    IF (d1 .le. d2   .and.   d1 .le. d3   .and.   d1 .le.d4 )   THEN 
        nearest_vertice_to_liquid_cell_cntr(1)= x_S10(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(2)= y_S10(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(3)= 1.   ! S10
        nearest_vertice_to_liquid_cell_cntr(4)= d1
        nearest_vertice_to_liquid_cell_cntr(5)= min(d2, d3, d4)


    ELSEIF (d2 .le. d1   .and.   d2 .le. d3   .and.   d2 .le.d4 )  THEN
        nearest_vertice_to_liquid_cell_cntr(1)= x_S01(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(2)= y_S01(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(3)= 2.   ! S01
        nearest_vertice_to_liquid_cell_cntr(4)= d2
        nearest_vertice_to_liquid_cell_cntr(5)= min(d1, d3, d4)

    ELSEIF (d3 .le. d1   .and.   d3 .le. d2   .and.   d3 .le.d4 ) THEN
        nearest_vertice_to_liquid_cell_cntr(1)= x_Sn10(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(2)= y_Sn10(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(3)= 3.   ! Sn10
        nearest_vertice_to_liquid_cell_cntr(4)= d3
        nearest_vertice_to_liquid_cell_cntr(5)= min(d1, d2, d4)


    ELSE
        nearest_vertice_to_liquid_cell_cntr(1)= x_S0n1(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(2)= y_S0n1(xx,yy)
        nearest_vertice_to_liquid_cell_cntr(3)= 4.   ! S10
        nearest_vertice_to_liquid_cell_cntr(4)= d4
        nearest_vertice_to_liquid_cell_cntr(5)= min(d1, d2, d3)

    ENDIF


END FUNCTION nearest_vertice_to_liquid_cell_cntr





REAL FUNCTION angle_between_liquid_cell_cntr_and_nearest_vertice(x,y,xx,yy,xS_near,yS_near)

USE m_dengrowth, ONLY: x_G, y_G

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: y, x, xx, yy
  REAL, INTENT(IN) :: xS_near, yS_near
  REAL, DIMENSION(2) ::  vector_G_S, vector_G_C
  REAL :: w_G_S, w_G_Cnt

    
vector_G_S= (/ xS_near-x_G(xx,yy), yS_near-y_G(xx,yy) /)
vector_G_C= (/ x*dx-x_G(xx,yy), y*dx-y_G(xx,yy) /)
w_G_S= sqrt ( (xS_near-x_G(xx,yy)) **2. +  (yS_near-y_G(xx,yy)) **2. )
w_G_Cnt= sqrt ( (x*dx-x_G(xx,yy)) **2. + (y*dx-y_G(xx,yy)) **2. )
angle_between_liquid_cell_cntr_and_nearest_vertice= acos (  dot_product(vector_G_S,vector_G_C) /  ( w_G_Cnt * w_G_S ) ) 

END FUNCTION angle_between_liquid_cell_cntr_and_nearest_vertice







END SUBROUTINE dendgrowth_CA



