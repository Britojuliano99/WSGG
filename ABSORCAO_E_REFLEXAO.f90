SUBROUTINE ABSORCAO_E_REFLEXAO(E_M,D_M,NANO_IN_VOLUME,CA_CINZA,CESP,CEXT,CESP_NANO,CABS_NANO,CEXT_NANO,R_N,FV,R_VC,FV_VC,&
&G,N_TE,CPART,CTEC,QA,QS,raio)

use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
INTEGER KK,ND22,N,NI,II,I,NZ,NR,NA,NCG,NUMR,NRZ1,NRZ2,NRR,NRG,ZV,RV,AV,ZV_E,RV_E,AV_E
INTEGER ZV_C,RV_C,AV_C,ZV_G,RV_G,AV_G,ZV_A,RV_A,AV_A
INTEGER NZ1,NZ1_LASER,NZ2,NUR,NG,NSZ1,NSZ1_LASER,NSZ2,NSR,NSG,PACOTE,DIR,MIST
INTEGER EXTINCTED_IN_NANO,BUNDLE_EXTINCTED,N_TE,I_TE,CPART,CTEC
INTEGER INT_TEMP

REAL*8 II1,TG,DIMZ,DIMR,DIM_ALPHA,DIMCG,DELZ,DELR,NP_Z,NP_R,NP,W,TREF,TE,FI,BE,GA,Z,R,AL
DOUBLE PRECISION Z_C,R_C,AL_C,Z_E,R_E,AL_E,Z_G,R_G,AL_G,ZV1,RV1,ZV_E1,RV_E1,ZV_C1,RV_C1,ZV_G1,RV_G1,K
DOUBLE PRECISION D,DN,DFI,DSEN,DCOS,DCOSI,DNIT1,DNIT2,DNI,DFII,FIIPT,FIRP
DOUBLE PRECISION DCOSRP2T,DRPT,DCOSRNT,DNRNT,DRNT,FIRPT,FIRPTI,TERP,BE_I,GA_I
DOUBLE PRECISION TERPI,DNRP,DFIRP,DSENI,DCOSIT,DCOSIT1,FI_I,TE_I
DOUBLE PRECISION X,Y,X_C,Y_C,C,DX,DY,DY1,DY2,DDY1,DDY2,DZ,DXY,K1,K2,K3,K4
DOUBLE PRECISION X1,X2,Y_C1,Y_C2,PI 
DOUBLE PRECISION PZ_T1,PZ_T2,PZ_T3,PZ_T4,PZ_T5,PZ_T6,PZ_T7,PZ_T8,PZ_T9 
DOUBLE PRECISION CE1,CE2,CE3,CEE,KC1,KC2,KC3,K_C,PW,PC,PWC,FPAG,K_ABS
DOUBLE PRECISION E0,R0,MU_T,FWHM,P_LASER
REAL*8 CABS_NANO,CEXT_NANO,CESP_NANO,R_N,FV,R_VC,FV_VC,NUMNANO_D
REAL*8 RAND_NANO,D_NANO,D_TEMP,RAND_NANO_D
REAL*8 L_C,TEI,G,D_TE,QA,QS,STE,RAIO_REAL,FV_REAL,TE_TEMP,TES,aaa

DIMENSION NANO_IN_VOLUME(NZ+2,NR+2,NA) 
DIMENSION E_M(NZ+2,NR+2,NA),D_M(NZ+2,NR+2,NA)
DIMENSION CESP(NZ+2,NR+2,NA),CEXT(NZ+2,NR+2,NA),CA_CINZA(NZ+2,NR+2,NA)
DIMENSION SRTEP(N_TE),SRTET(N_TE),SRTEP_MIE(N_TE)

character(50) raio


COMMON II1,TG,DIMZ,DIMR,DIM_ALPHA,DIMCG,DELZ,DELR,NP_Z,NP_R,NP,W,TREF,TE,FI,BE,GA,Z,R,AL
COMMON Z_C,R_C,AL_C,Z_E,R_E,AL_E,Z_G,R_G,AL_G,ZV1,RV1,ZV_E1,RV_E1,ZV_C1,RV_C1,ZV_G1,RV_G1,K
COMMON D,DN,DFI,DSEN,DCOS,DCOSI,DNIT1,DNIT2,DNI,DFII,FIIPT,FIRP
COMMON DCOSRP2T,DRPT,DCOSRNT,DNRNT,DRNT,FIRPT,FIRPTI,TERP,BE_I,GA_I
COMMON TERPI,DNRP,DFIRP,DSENI,DCOSIT,DCOSIT1,FI_I,TE_I
COMMON X,Y,X_C,Y_C,C,DX,DY,DY1,DY2,DDY1,DDY2,DZ,DXY,K1,K2,K3,K4
COMMON X1,X2,Y_C1,Y_C2,PI 
COMMON PZ_T1,PZ_T2,PZ_T3,PZ_T4,PZ_T5,PZ_T6,PZ_T7,PZ_T8,PZ_T9 
COMMON CE1,CE2,CE3,CEE,KC1,KC2,KC3,K_C,PW,PC,PWC,FPAG,K_ABS
COMMON E0,R0,MU_T,FWHM,P_LASER

COMMON KK,ND22,N,NI,II,I,NZ,NR,NA,NCG,NUMR,NRZ1,NRZ2,NRR,NRG,ZV,RV,AV,ZV_E,RV_E,AV_E
COMMON ZV_C,RV_C,AV_C,ZV_G,RV_G,AV_G,ZV_A,RV_A,AV_A
COMMON NZ1,NZ1_LASER,NZ2,NUR,NG,NSZ1,NSZ1_LASER,NSZ2,NSR,NSG,PACOTE,DIR,MIST


IF (CPART == 1) THEN 

   D_TE = PI/(N_TE-1)

   TES = 0.D0

   SRTEP(1) =  SIN(TES)*D_TE/2.d0 

   DO I_TE = 2,N_TE 

      TES = TES + D_TE

      SRTEP(I_TE) = SRTEP(I_TE-1) + SIN(TES)*D_TE/2.d0 

   ENDDO

ENDIF


IF (CPART == 2) THEN 

   D_TE = PI/(N_TE-1)

   TES = 0.D0

   SRTEP(1) = ((1.D0 - G**2.D0)/((1.D0 + G**2.D0 - 2.D0*G*COS(TES))**(1.5D0)))*SIN(TES)*D_TE/2.D0

   DO I_TE = 2,N_TE 

      TES = TES + D_TE

      SRTEP(I_TE) = SRTEP(I_TE-1) + ((1.D0 - G**2.D0)/((1.D0 + G**2.D0 - 2.D0*G*COS(TES))**(1.5D0)))*SIN(TES)*D_TE/2.D0

   ENDDO

ENDIF

IF (CPART == 3) THEN 

   D_TE = PI/(N_TE-1)

   TES = 0.D0

   SRTEP(1) = 0.75D0*(1.D0 + COS(TES)**2.D0)*SIN(TES)*D_TE/2.D0 

   DO I_TE = 2,N_TE 

      TES = TES + D_TE

      SRTEP(I_TE) = SRTEP(I_TE-1) + 0.75D0*(1.D0 + COS(TES)**2.D0)*SIN(TES)*D_TE/2.D0 

   ENDDO

ENDIF

IF (CPART == 4) THEN 

   OPEN(333,FILE='phase_function/phi_Ra615_Rb75_L6328.txt')

      DO I_TE=1,N_TE

         READ(333,*)aaa,SRTEP_MIE(I_TE)

      ENDDO

   CLOSE(333)

   D_TE = PI/(N_TE-1)

   TES = 0.D0

   SRTEP(1) = SIN(TES)*SRTEP_MIE(1)*D_TE/2.d0

   DO I_TE = 2,N_TE 

      TES = TES + D_TE

      SRTEP(I_TE) = SRTEP(I_TE-1) + SRTEP_MIE(I_TE)*SIN(TES)*D_TE/2.d0  

   ENDDO

ENDIF

!OPEN(111,FILE='/home/labsistermicos/Documentos/Andre/MonteCarlo/PTT_WSGG_GC_3D_SD/phi_plot.m')

!write(111,*)'phi_temp=['
!      DO I_TE=2,N_TE

!         write(111,*)SRTEP(I_TE)

!      ENDDO
!write(111,*)'];'
!CLOSE(111)
!pause

IF (CTEC == 1) THEN 

   D_TE = PI/(N_TE-1)

   TES = 0.D0

   SRTET(1) =  SIN(TES)*D_TE/2.d0 

   DO I_TE = 2,N_TE 

      TES = TES + D_TE

      SRTET(I_TE) = SRTET(I_TE-1) + SIN(TES)*D_TE/2.d0 

   ENDDO

ENDIF


IF (CTEC == 2) THEN 

   D_TE = PI/(N_TE-1)

   TES = 0.D0

   SRTET(1) = ((1.D0 - G**2.D0)/((1.D0 + G**2.D0 - 2.D0*G*COS(TES))**(1.5D0)))*SIN(TES)*D_TE/2.D0

   DO I_TE = 2,N_TE 

      TES = TES + D_TE

      SRTET(I_TE) = SRTET(I_TE-1) + ((1.D0 - G**2.D0)/((1.D0 + G**2.D0 - 2.D0*G*COS(TES))**(1.5D0)))*SIN(TES)*D_TE/2.D0

   ENDDO

ENDIF

IF (CTEC == 3) THEN 

   D_TE = PI/(N_TE-1)

   TES = 0.D0

   SRTET(1) = 0.75D0*(1.D0 + COS(TES)**2.D0)*SIN(TES)*D_TE/2.D0

   DO I_TE = 2,N_TE 

      TES = TES + D_TE

      SRTET(1) = 0.75D0*(1.D0 + COS(TES)**2.D0)*SIN(TES)*D_TE/2.D0  

   ENDDO

ENDIF

 

NUMNANO_D = (3.D0*FV/(4.d0*PI*R_N**3.d0))**(1.d0/3.d0)  !NUMBER OF NANOPARTICLES PER METER  


INT_TEMP = INT(NR*R_E/DIMR)

IF (R_E - INT_TEMP*DIMR/NR .GT. 0.D0 ) THEN 

   RV_E = INT_TEMP + 2

   ELSE 

   RV_E = INT_TEMP + 1

   IF (R_E .EQ. DIMR) THEN 

      RV_E = NR + 2
    
   ENDIF

ENDIF

INT_TEMP = INT(NZ*Z_E/DIMZ)

IF (Z_E - INT_TEMP*DIMZ/NZ .GT. 0.D0 ) THEN

   ZV_E = INT_TEMP + 2

   ELSE

   ZV_E = INT_TEMP + 1

   IF (Z_E .EQ. DIMZ) THEN 

      ZV_E = NZ + 2
    
   ENDIF

ENDIF

 
INT_TEMP =  INT(NA*AL_E/DIM_ALPHA)

IF (AL_E - INT_TEMP*DIM_ALPHA/NA .GT. 0.D0 ) THEN 

   AV_E =  INT_TEMP + 1

   ELSE 

   AV_E =  INT_TEMP 

ENDIF



DO WHILE (PACOTE==1) 

   INT_TEMP = INT(NR*R_C/DIMR)

   IF (R_C - INT_TEMP*DIMR/NR .GT. 0.D0 ) THEN 

      RV_C = INT_TEMP + 2

      ELSE 

      RV_C = INT_TEMP + 1

      IF (R_C .EQ. DIMR) THEN 

         RV_C = NR + 2
    
      ENDIF

   ENDIF

   INT_TEMP = INT(NZ*Z_C/DIMZ)

   IF (Z_C - INT_TEMP*DIMZ/NZ .GT. 0.D0 ) THEN 

      ZV_C = INT_TEMP + 2

      ELSE

      ZV_C = INT_TEMP + 1

      IF (Z_C .EQ. DIMZ) THEN 

         ZV_C = NZ + 2
    
      ENDIF

   ENDIF

 
   INT_TEMP =  INT(NA*AL_C/DIM_ALPHA)

   IF (AL_C - INT_TEMP*DIM_ALPHA/NA .GT. 0.D0 ) THEN 

      AV_C =  INT_TEMP + 1

      ELSE 

      AV_C =  INT_TEMP 

   ENDIF

 !     IF (ZV_C < 0.D0 .OR. RV_C < 0.D0) THEN

 !     write(*,*)'FIRP,TERP',FIRP,TERP,Z_C,R_C,DZ

 !     write(*,*)'2'X,Y,C,X1,X2,X_C,DXY,TERP,TAN(TERP)

 !     write(*,*)'3',DX,DY

 !     ENDIF

   D=DIMCG/NCG

   DAUXM=((X_C-X)**2.D0+(Y_C-Y)**2.D0+(Z_C-Z)**2.D0)**.5D0
   
  
   
   ND=NINT(DAUXM/D)

   
   IF (ND==0) THEN
      ND=1
   ENDIF

   ND=2*ND
   D=DAUXM/ND 

   ND2=0
   DZG=D*SIN(TERP) 
   DXYG=D*COS(TERP)
   DXG=DXYG*COS(FIRP)
   DYG=DXYG*SIN(FIRP)



   X_G=X
   Y_G=Y
   Z_G=Z

   PROB=2.D0
   EXTINCTED_IN_NANO = 0 !EXTINCTED_IN_NANO = 1 MEANS THAT THE BUNDLE WAS EXTINCTED IN THE NANOPARTICLE (SCATTERED OR ABSORBED).   

   BUNDLE_EXTINCTED = 0 ! BUNDLE_EXTINCTED = 1 MEANS THAT THE BUNDLE WAS EXTINCTED (SCATTERED OR ABSORBED). 


   DO WHILE (BUNDLE_EXTINCTED == 0 .AND. ND2 < ND)

      D=DAUXM/ND 

      D_TEMP = D

      ND2=ND2+1


      X_G=X_G+DXG
      Y_G=Y_G+DYG
      Z_G=Z_G+DZG
      R_G=((X_G-DIMR)**2.D0+(Y_G-DIMR)**2.D0)**.5D0
      !AL_G=ATAN((DIMR-Y_G)/(DIMR-X_G))


       IF (Y_G .GT. DIMR) THEN
         
         IF (X_G .GT. DIMR) THEN

            AL_G=ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

         IF (X_G .LT. DIMR) THEN

            AL_G= PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

      ENDIF


      IF (Y_G .LT. DIMR) THEN
         
         IF (X_G .GT. DIMR) THEN

            AL_G= 2.D0*PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

         IF (X_G .LT. DIMR) THEN

            AL_G= PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

      ENDIF

      IF (X_G .EQ. DIMR) THEN 

         IF (Y_G .GE. DIMR) THEN 

            AL_G = PI/2.D0

         ENDIF

         IF (Y_G .LT. DIMR) THEN 

            AL_G = 3.D0*PI/2.D0

         ENDIF

      ENDIF

      IF (Y_G .EQ. DIMR) THEN 

         IF (X_G .GE. DIMR) THEN 

            AL_G = 0.D0

         ENDIF

         IF (X_G .LT. DIMR) THEN 

            AL_G = PI

         ENDIF

      ENDIF


      !   AL_G = PI + AL_G
      !   ELSE
      !   IF (Y_G>DIMR ) THEN

      !      AL_G = 2.D0*PI + AL_G
      !   ENDIF
      !ENDIF      



      !IF (X_G>DIMR) THEN
      !   AL_G = PI + AL_G
      !   ELSE
      !   IF (Y_G>DIMR ) THEN
      !      AL_G = 2.D0*PI + AL_G
      !   ENDIF
      !ENDIF
     

      INT_TEMP = INT(NR*R_G/DIMR)

      IF (R_G - INT_TEMP*DIMR/NR .GT. 0.D0 ) THEN 

         RV_G = INT_TEMP + 2

         ELSE 

         RV_G = INT_TEMP + 1

      ENDIF

      INT_TEMP = INT(NZ*Z_G/DIMZ)

      IF (Z_G - INT_TEMP*DIMZ/NZ .GT. 0.D0 ) THEN 

         ZV_G = INT_TEMP + 2

         ELSE

         ZV_G = INT_TEMP + 1

      ENDIF

!AL_G = 2.d0*2.d0*PI/3.d0 - 1.d-8

 
      INT_TEMP =  INT(NA*AL_G/DIM_ALPHA)

      IF (AL_G - INT_TEMP*DIM_ALPHA/NA .GT. 0.D0 ) THEN 

         AV_G =  INT_TEMP + 1

         ELSE 

         AV_G =  INT_TEMP 

      ENDIF

      IF (Z_G .LE. 0.D0) THEN
      
         ZV_G = 2

      ENDIF 

      IF (Z_G .GE. DIMZ) THEN
      
         ZV_G = NZ + 1

      ENDIF 

      IF (R_G .GE. DIMR) THEN
      
         RV_G = NR + 1

      ENDIF 

!write(*,*)'AL_G,AV_G',AL_G,AV_G  

!read(stdin,*)


      IF (ZV_G > 0.D0 .AND. RV_C > 0.D0) THEN

      IF (NANO_IN_VOLUME(ZV_G,RV_G,AV_G) == 1) THEN ! THERE ARE NANOPATICLES IN THE VOLUME.

         D_TEMP = D 

         !!write(*,*)'D_TEMP',D_TEMP

         DO WHILE (D_TEMP > 0.D0) 
  
            CALL RANDOM_NUMBER(RAND_NANO)
       
            D_NANO = (1.D0/NUMNANO_D)*LOG(1.D0-RAND_NANO)/LOG(1.D0-(NUMNANO_D**2.d0)*(PI*R_VC**2.d0)) ! LEGTH TRAVELED BY THE BUNDLE BEFORE REACHING A NANOPARTICLE.
            
            
            IF (D_NANO < D_TEMP) THEN  ! THERE IS A NANOPARTICLE IN D. 

               CALL RANDOM_NUMBER(R_PG)
 
               PG=1.D0-EXP(-D_NANO*CEXT(ZV_G,RV_G,AV_G)) 

               IF (R_PG < PG) THEN ! BUNDLE IS EXTINCTED BEFORE REACHING THE NANOPARTICLE.

                  EXTINCTED_IN_NANO = 0

                  BUNDLE_EXTINCTED = 1 
     
                  L_ESP = D/2.D0

                  D_TEMP = 0.D0 

               ELSE                 ! BUNDLE REACHES THE NANOPARTICLE.

               
                  CALL RANDOM_NUMBER(R_PG)

                  PG=  1.d0 !AQUI 1.D0-EXP(-L_C*CEXT_NANO)
                       
                 
                  IF (R_PG < PG) THEN  ! BUNDLE IS EXTINCTED IN THE NANOPARTICLE. 
 
                     EXTINCTED_IN_NANO = 1

                     BUNDLE_EXTINCTED = 1 
     
                     L_ESP = D/2.D0

                     D_TEMP = 0.D0 

                  ELSE

                     IF (D_NANO < D_TEMP) THEN 

                        D_TEMP = D_TEMP - D_NANO 

                     ELSE 

                        D_TEMP = 0.D0 

                     ENDIF 

                  ENDIF

               ENDIF

            ELSE

               CALL RANDOM_NUMBER(R_PG)

               PG=1.D0-EXP(-D_TEMP*CEXT(ZV_G,RV_G,AV_G)) 

               IF (R_PG < PG) THEN 
           
                  EXTINCTED_IN_NANO = 0

                  BUNDLE_EXTINCTED = 1 
     
                  L_ESP = D/2.D0

                  D_TEMP = 0.D0 
               
               ENDIF  

                  D_TEMP = 0.D0 

            ENDIF

         ENDDO

      ENDIF

      IF (NANO_IN_VOLUME(ZV_G,RV_G,AV_G) == 0) THEN ! THERE ARE NOT NANOPATICLES IN THE VOLUME.

          CALL RANDOM_NUMBER(R_PG)

         PG=1.D0-EXP(-D*CEXT(ZV_G,RV_G,AV_G))   

         IF (R_PG < PG) THEN 

            EXTINCTED_IN_NANO = 0

            BUNDLE_EXTINCTED = 1 
     
            L_ESP = D/2.D0
                
         ENDIF

      ENDIF

      ENDIF

   ENDDO

   IF (ZV_G > 0.D0 .AND. RV_G > 0.D0) THEN


   IF (BUNDLE_EXTINCTED == 1) THEN !THE BUNDLE IS EITHER ABSORBED OR SCATTERD AWAY FROM THE PROPAGATION PATH.

      IF (EXTINCTED_IN_NANO == 1) THEN

         PABS = QA/(QA+QS)

         CALL RANDOM_NUMBER(R_ABS) 


         IF (R_ABS <= PABS) THEN ! BUNDLE IS ABSORVIDO. 

   	    ZV_A=ZV_G

	    RV_A=RV_G

            AV_A = AV_G

            PACOTE=0
    
         !!write(*,*)'BUNDLE IS ABSORVIDO'

         ELSE                    ! BUNDLE IS SCATTERED.


         X_G=X_G-DXG
         Y_G=Y_G-DYG
         Z_G=Z_G-DZG

      R_G=((X_G-DIMR)**2.D0+(Y_G-DIMR)**2.D0)**.5D0
      

      IF (Y_G .GT. DIMR) THEN
         
         IF (X_G .GT. DIMR) THEN

            AL_G=ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

         IF (X_G .LT. DIMR) THEN

            AL_G= PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

      ENDIF


      IF (Y_G .LT. DIMR) THEN
         
         IF (X_G .GT. DIMR) THEN

            AL_G= 2.D0*PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

         IF (X_G .LT. DIMR) THEN

            AL_G= PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

      ENDIF

      IF (X_G .EQ. DIMR) THEN 

         IF (Y_G .GE. DIMR) THEN 

            AL_G = PI/2.D0

         ENDIF

         IF (Y_G .LT. DIMR) THEN 

            AL_G = 3.D0*PI/2.D0

         ENDIF

      ENDIF

      IF (Y_G .EQ. DIMR) THEN 

         IF (X_G .GE. DIMR) THEN 

            AL_G = 0.D0

         ENDIF

         IF (X_G .LT. DIMR) THEN 

            AL_G = PI

         ENDIF

      ENDIF





      !AL_G=ATAN((DIMR-Y_G)/(DIMR-X_G))
      !IF (X_G>DIMR) THEN
      !   AL_G = PI + AL_G
      !   ELSE
      !   IF (Y_G>DIMR ) THEN
      !      AL_G = 2.D0*PI + AL_G
      !   ENDIF
      !ENDIF
  
      INT_TEMP = INT(NR*R_G/DIMR)

      IF (R_G - INT_TEMP*DIMR/NR .GT. 0.D0 ) THEN 

         RV_G = INT_TEMP + 2

         ELSE 

         RV_G = INT_TEMP + 1

      ENDIF

      INT_TEMP = INT(NZ*Z_G/DIMZ)

      IF (Z_G - INT_TEMP*DIMZ/NZ .GT. 0.D0 ) THEN 

         ZV_G = INT_TEMP + 2

         ELSE

         ZV_G = INT_TEMP + 1

      ENDIF


      INT_TEMP =  INT(NA*AL_G/DIM_ALPHA)

      IF (AL_G - INT_TEMP*DIM_ALPHA/NA .GT. 0.D0 ) THEN 

         AV_G =  INT_TEMP + 1

         ELSE 

         AV_G =  INT_TEMP 

      ENDIF

      IF (Z_G .LE. 0.D0) THEN
      
         ZV_G = 2

      ENDIF 

      IF (Z_G .GE. DIMZ) THEN
      
         ZV_G = NZ + 1

      ENDIF 

      IF (R_G .GE. DIMR) THEN
      
         RV_G = NR + 1

      ENDIF 



      IF (ZV_G > 0.D0 .AND. RV_C > 0.D0) THEN

            CALL RANDOM_NUMBER(R_TE)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Angle of scattering by the nanoparticle
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            I_TE = 1            
            TE_TEMP = TE
            TE = 0.D0
            DO WHILE (SRTEP(I_TE) < R_TE .AND. I_TE < N_TE)
               TE = TE + D_TE
               I_TE = I_TE + 1
            ENDDO               
            TE = TE_TEMP + TE         
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
            CALL RANDOM_NUMBER(R_FI)
            FI=2.D0*PI*R_FI
            IF (R_FI<.5D0) THEN
               NFI=NFI+1
            ENDIF
            BE=0.D0
            GA=0.D0
            Z = Z_G 
            R = R_G 
            AL = AL_G  
           
            CALL CHEGADA      
    
         ENDIF

         ENDIF
          
        
      ELSE

        PABS= (1.D0 - EXP(-D*CA_CINZA(ZV_G,RV_G,AV_G)))/(1.D0-EXP(-D*CEXT(ZV_G,RV_G,AV_G)))

        CALL RANDOM_NUMBER(R_ABS) 

        IF (R_ABS <= PABS) THEN ! BUNDLE IS ABSORVIDO. 

     
   	    ZV_A=ZV_G

	    RV_A=RV_G

            AV_A = AV_G

            PACOTE=0

         ELSE                    ! BUNDLE IS SCATTERED.

            X_G=X_G-DXG
            Y_G=Y_G-DYG
            Z_G=Z_G-DZG


       R_G=((X_G-DIMR)**2.D0+(Y_G-DIMR)**2.D0)**.5D0




      IF (Y_G .GT. DIMR) THEN
         
         IF (X_G .GT. DIMR) THEN

            AL_G=ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

         IF (X_G .LT. DIMR) THEN

            AL_G= PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

      ENDIF


      IF (Y_G .LT. DIMR) THEN
         
         IF (X_G .GT. DIMR) THEN

            AL_G= 2.D0*PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

         IF (X_G .LT. DIMR) THEN

            AL_G= PI + ATAN((DIMR-Y_G)/(DIMR-X_G))

         ENDIF

      ENDIF

      IF (X_G .EQ. DIMR) THEN 

         IF (Y_G .GE. DIMR) THEN 

            AL_G = PI/2.D0

         ENDIF

         IF (Y_G .LT. DIMR) THEN 

            AL_G = 3.D0*PI/2.D0

         ENDIF

      ENDIF

      IF (Y_G .EQ. DIMR) THEN 

         IF (X_G .GE. DIMR) THEN 

            AL_G = 0.D0

         ENDIF

         IF (X_G .LT. DIMR) THEN 

            AL_G = PI

         ENDIF

      ENDIF


      !AL_G=ATAN((DIMR-Y_G)/(DIMR-X_G))
      !IF (X_G>DIMR) THEN
      !   AL_G = PI + AL_G
      !   ELSE
      !   IF (Y_G>DIMR ) THEN
      !      AL_G = 2.D0*PI + AL_G
      !   ENDIF
      !ENDIF
     

      INT_TEMP = INT(NR*R_G/DIMR) + 1

      IF (R_G - INT_TEMP*DIMR/NR .LE. 0.D0 ) THEN 

         RV_G = INT_TEMP + 1

         ELSE 

         RV_G = INT_TEMP 

      ENDIF


      INT_TEMP = INT(NZ*Z_G/DIMZ) + 1

      IF (Z_G - INT_TEMP*DIMZ/NZ .LE. 0.D0 ) THEN 

         ZV_G = INT_TEMP + 1

         ELSE

         ZV_G = INT_TEMP 

      ENDIF


      INT_TEMP =  INT(NA*AL_G/DIM_ALPHA)

      IF (AL_G - INT_TEMP*DIM_ALPHA/NA .GT. 0.D0 ) THEN 

         AV_G =  INT_TEMP + 1

         ELSE 

         AV_G =  INT_TEMP 

      ENDIF



   ! INT_TEMP =  INT(NA*AL_G/DIM_ALPHA)

   !IF (AL_G - INT_TEMP*DIM_ALPHA/NA .LE. 0.D0 ) THEN 

    !  AV_G =  INT_TEMP

    !  ELSE 

    !  AV_G =  INT_TEMP + 1

   !ENDIF

   IF (ZV_G > 0.D0 .AND. RV_C > 0.D0) THEN

            CALL RANDOM_NUMBER(R_TE)

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !Angle of scattering by the tissue
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            I_TE = 1               
            TE_TEMP = TE
            TE = 0.D0
            DO WHILE (SRTET(I_TE) < R_TE .AND. I_TE < N_TE) 
               TE = TE + D_TE
               I_TE = I_TE + 1
            ENDDO                      
            TE = TE_TEMP + TE
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            CALL RANDOM_NUMBER(R_FI)
            FI=2.D0*PI*R_FI
            IF (R_FI<.5D0) THEN
               NFI=NFI+1
            ENDIF
            BE=0.D0
            GA=0.D0
            Z = Z_G 
            R = R_G 
            AL = AL_G  
           
            CALL CHEGADA          

  ENDIF
         ENDIF

      ENDIF

   ENDIF

   IF (ZV_G > 0.D0 .AND. RV_C > 0.D0) THEN

   IF (BUNDLE_EXTINCTED == 0) THEN ! BUNDLE REACHS A SURFACE.


      CALL RANDOM_NUMBER(R_ALFA)
      
      ALFA=E_M(ZV_C,RV_C,AV_C)

      IF (ALFA >= R_ALFA) THEN ! BUNDLE IS ABSORBED BY A SURFACE.

         ZV_A=ZV_C

	 RV_A=RV_C

         AV_A=AV_C

         PACOTE=0

 
      ENDIF

      IF (ALFA < R_ALFA) THEN ! BUNDLE IS REFLECTED BY A SURFACE. 

         IF (D_M(ZV_C,RV_C,AV_C) > 0) THEN ! REFLECTION IS 100% SPECULAR. 

            TE = TE_I   

	    FI=FI_I+PI

            Z=Z_C

            R=R_C

            AL=AL_C

            BE=BE_I

            GA=GA_I

            CALL CHEGADA

	 ELSE                         ! REFLECTION IS 100% DIFFUSE. 
     
	    CALL RANDOM_NUMBER(R_TE)

	    CALL RANDOM_NUMBER(R_FI)

            TE=PI/2.D0-ASIN(R_TE**.5D0)

            FI=2.D0*PI*R_FI

            Z=Z_C

            R=R_C

            AL=AL_C

            BE=BE_I

            GA=GA_I

            CALL CHEGADA

	 ENDIF    
      
      ENDIF
   
    ENDIF 
  
    ENDIF 

    ENDIF

    IF (ZV_G .LE. 0.d0 .OR. RV_G .LE. 0.d0) THEN

       PACOTE = 0

    ENDIF

ENDDO

END SUBROUTINE ABSORCAO_E_REFLEXAO

