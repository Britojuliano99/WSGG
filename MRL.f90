SUBROUTINE MRL(Z_E,R_E,AL_E,Z,R,AL,Z_LASER,R_LASER,AL_LASER,BE_LASER,GA_LASER,TE,FI)

! MRL: Mundança de Referencial do Laser.

REAL*8 Z_E,R_E,AL_E,Z,R,AL,Z_LASER,R_LASER,AL_LASER,BE_LASER,GA_LASER,TE,FI
REAL*8 X1,X2,X3,Y1,Y2,Y3,Z1,Z2,Z3,PI,XN,YN,ZN,H,H2



PI = acos(-1.d0)


X1 = R_E*cos(AL_E)

Y1 = R_E*sin(AL_E)


Y2 = Y1*cos(BE_LASER)

X2 = X1*cos(GA_LASER)

Z2 = Y1*sin(BE_LASER) 


X3 = R_LASER*cos(AL_LASER) + X2 

Y3 = R_LASER*sin(AL_LASER) + Y2

Z3 = Z_LASER + Z2

Z = Z3

R = (X3**2.D0 + Y3**2.D0)**.5D0


!write(*,*)'Y3',Y3

!pause

IF (Y3 >= 0.D0) THEN

   IF (X3 >= 0.d0) THEN 

      ! Primeiro quadrante
      
      AL = atan(Y3/X3) 

      ELSE

      ! Segundo quadrante 

      AL = PI + atan(Y3/X3) 

   ENDIF

   ELSE

   IF (X3 < 0.D0 ) THEN 

      ! Terceiro quadrante

      AL = PI + atan(Y3/X3)   

      ELSE 

      ! Quarto quadrante

      AL = 2.D0*PI + atan(Y3/X3) 

   ENDIF

ENDIF 


IF (X3 .EQ. 0.d0) THEN 

   IF (Y3 .GE. 0.D0) THEN

              
      AL = PI/2.D0 

   ENDIF

   IF (Y3 .LT. 0.D0) THEN

      AL = 3.D0*PI/2.D0 

   ENDIF

ENDIF


IF (Y3 .EQ. 0.d0) THEN 

   IF (X3 .GE. 0.D0) THEN

              
      AL = 0.D0 

   ENDIF

   IF (X3 .LT. 0.D0) THEN

      AL = PI 

   ENDIF

ENDIF

  




!IF (GA_LASER .NE. 0.D0 .OR. BE_LASER .NE. 0.D0) THEN

XN = sin(GA_LASER)

YN = sin(BE_LASER)

ZN = cos(BE_LASER)

H = (XN**2.D0 + YN**2.D0)**.5d0

H2 = (XN**2.D0 + YN**2.D0 + ZN**2.D0)**.5d0

TE = acos(H/H2)*(ZN/abs(ZN))


   IF (YN > 0.D0) THEN

      IF (XN > 0.D0) THEN 
         
         !Primeiro quadrante

         FI = atan(YN/XN) 

      ENDIF

      IF (XN < 0.D0) then 

         !Segundo quadrante 

         FI = PI + atan(YN/XN) 

      ENDIF

   ENDIF

   IF (YN < 0.D0) THEN

      IF (XN < 0.D0) THEN 

         ! Terceiro quadrante

         FI = PI + atan(YN/XN)   

      ENDIF

      IF (XN < 0.D0) THEN 

         ! Quarto quadrante

         FI = 2.D0*PI + atan(YN/XN) 

      ENDIF

   ENDIF
 
!ENDIF

  
END SUBROUTINE MRL
