PROGRAM CILINDRO_PPT

INTEGER KK,ND22,N,NI,II,I,NZ,NR,NA,NCG,NUMR,NRZ1,NRZ2,NRR,NRG,ZV,RV,AV,ZV_E,RV_E,AV_E
INTEGER ZV_C,RV_C,AV_C,ZV_G,RV_G,AV_G,ZV_A,RV_A,AV_A
INTEGER NZ1,NZ1_LASER,NZ2,NUR,NG,NSZ1,NSZ1_LASER,NSZ2,NSR,NSG,PACOTE,DIR,MIST
INTEGER F_MULT,N_TE,CROSS_SECTION
INTEGER NR_LASER1,seed,CPART,CTEC,Domb,N_DESV,I_DESV,N_LASER1

REAL*8 II1,TG,DIMZ,DIMR,DIMCG,DELZ,DELR,NP_Z,NP_R,NP,W,TREF,TE,FI,BE,GA,Z,R,AL
DOUBLE PRECISION Z_C,R_C,AL_C,Z_E,R_E,AL_E,Z_G,R_G,AL_G,ZV1,RV1,ZV_E1,RV_E1,ZV_C1,RV_C1,ZV_G1,RV_G1,K
DOUBLE PRECISION D,DN,DFI,DSEN,DCOS,DCOSI,DNIT1,DNIT2,DNI,DFII,FIIPT,FIRP
DOUBLE PRECISION DCOSRP2T,DRPT,DCOSRNT,DNRNT,DRNT,FIRPT,FIRPTI,TERP,BE_I,GA_I
DOUBLE PRECISION TERPI,DNRP,DFIRP,DSENI,DCOSIT,DCOSIT1,FI_I,TE_I
DOUBLE PRECISION X,Y,X_C,Y_C,C,DX,DY,DY1,DY2,DDY1,DDY2,DZ,DXY,K1,K2,K3,K4
DOUBLE PRECISION X1,X2,Y_C1,Y_C2,PI

DOUBLE PRECISION PZ_T1,PZ_T2,PZ_T3,PZ_T4,PZ_T5,PZ_T6,PZ_T7,PZ_T8,PZ_T9
DOUBLE PRECISION CE1,CE2,CE3,CEE,KC1,KC2,KC3,K_C,PW,PC,PWC,FPAG,K_ABS!,CESP,CEXT,CA_CINZA
DOUBLE PRECISION E0,R0,MU_T,FWHM,P_LASER
REAL*8 QA,QS,QEXT,NN,NOVO1,NOVO2,RN2,ACROSS_A,ACROSS_S,ACROSS_E

REAL*8 CABS_NANO, CEXT_NANO,CESP_NANO,R_N,FV,R_VC,FV_VC
REAL*8 DIMR_LASER1,DIM_ALPHA,G,DESREDP,DESREDT
REAL*8 DESV,SDESV,MEAN,CONFIDENCE
REAL*8 Z_LASER1,R_LASER1,A_LASER1,Q_LASER1,AL_LASER1,GA_LASER1,BE_LASER1


character(50) pindex,tindex,raio

INTEGER,DIMENSION(:,:,:), ALLOCATABLE:: N_G_M,N_R_M,N_Z1_M,N_Z2_M,N_M_LASER
INTEGER,DIMENSION(:,:,:), ALLOCATABLE:: NANO_IN_VOLUME
REAL,DIMENSION(:), ALLOCATABLE:: DDZ,P_Z,P_A,SAMPLE
REAL,DIMENSION(:), ALLOCATABLE:: DDR,R_M,P_R
REAL,DIMENSION(:), ALLOCATABLE:: taxa_total 
REAL,DIMENSION(:,:), ALLOCATABLE:: QWX1,QWX2,QWX3
REAL,DIMENSION(:,:), ALLOCATABLE:: TAXA1,TAXA2,TAXA3
REAL,DIMENSION(:,:,:), ALLOCATABLE:: E_M,D_M
REAL,DIMENSION(:,:,:), ALLOCATABLE:: A_M,V_M,EMIT_M,EAB2,T_M
REAL,DIMENSION(:,:,:), ALLOCATABLE:: QBG
REAL,DIMENSION(:,:,:), ALLOCATABLE:: Q_LASER
REAL,DIMENSION(:,:,:), ALLOCATABLE:: CESP,CEXT,CA_CINZA,QG


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




WRITE(*,*)'WSGG_2D=CERTO'

WRITE(*,*)'CILINDRO_WSGG_2D RODANDO'
WRITE(*,*)''

CALL RANDOM_SEED

!***********************************************************************
!DADOS QUE PODEM SER MODIFICADOS
!***********************************************************************
!WRITE(*,*)'WSGG',PACOTE
!PAUSE   

PI = ACOS(-1.D0)

CONT1 = 0
CONT2 = 0 

N_DESV = 1


!DIMENSOES E NUMERO DE ZONAS
DIMZ =  1.D0 !5.0d0 !1.0D-3   !.2D0
DIMR =  .2D0 !5.D-3 !20.0D-9   !100.D0
DIM_ALPHA = 2.D0*PI
DIMR_LASER1 = .2D0 !DIMR !RAIO DO LASER 

Z_LASER1 =  0.d0
R_LASER1 =  0.D0
AL_LASER1 = 0.D0
BE_LASER1 = 0.d0
GA_LASER1 = 0.D0

F_MULT = 1
NZ = 30 
NR = 10 
NA = 1

A_LASER1 = PI*DIMR_LASER1**2.D0

!TAXA DE EMISSAO DO LASER [W] E NUMERO DE PACOTES DE ENERGIA EMITIDOS PELO LASER
Q_LASER1 = 0.d0 !2.d4*A_LASER1
N_LASER1 = 0


!FATOR DE PRECISAO NO CALCULO DA ABSORCAO NO CAMINHO DE GAS
!QUANTO MAIOR O VALOR DE FAPG MAIOR A PRECISAO
!OBS: GERALMENTE FAPG=1 GARENTE BOA PRECISAO
FPAG=1.D0 !CHANGES THE VALUE OF D (A PATH LENGTH STEP SIZE)

ALLOCATE(DDZ(NZ+2),P_Z(NZ+2))
ALLOCATE(DDR(NR+2),R_M(NR+2),P_R(NR+2),P_A(NA))
ALLOCATE(N_G_M(NZ+2,NR+2,NA),N_R_M(NZ+2,NR+2,NA),N_Z1_M(NZ+2,NR+2,NA))
ALLOCATE(N_Z2_M(NZ+2,NR+2,NA),N_M_LASER(NZ+2,NR+2,NA))
ALLOCATE(E_M(NZ+2,NR+2,NA),D_M(NZ+2,NR+2,NA))
ALLOCATE(A_M(NZ+2,NR+2,NA),V_M(NZ+2,NR+2,NA),EMIT_M(NZ+2,NR+2,NA))
ALLOCATE(EAB2(NZ+2,NR+2,NA),T_M(NZ+2,NR+2,NA))
ALLOCATE(QBG(NZ+2,NR+2,NA),QG(NZ,NR,NA))
ALLOCATE(Q_LASER(NZ+2,NR+2,NA))
ALLOCATE(CESP(NZ+2,NR+2,NA),CEXT(NZ+2,NR+2,NA),CA_CINZA(NZ+2,NR+2,NA))
!ALLOCATE(E_EMAB(NZ+2,NR+2,NA,NZ+2,NR+2,NA))
ALLOCATE(NANO_IN_VOLUME(NZ+2,NR+2,NA))
ALLOCATE(SAMPLE(N_DESV))
ALLOCATE(taxa_total(N_DESV))
ALLOCATE(QWX1(NR,NA),QWX2(NR,NA),QWX3(NZ,NA))
ALLOCATE(TAXA1(NR,NA),TAXA2(NR,NA),TAXA3(NZ,NA))


!EMISSIVIDADES E CARACTERISTICAS DE REFLEXAO (DIFUSA OU ESPECULAR)
!DAS SUPERFICIES
DO RV=1,NR+2
   DO AV = 1,NA
      E_M(1,RV,AV)= .5D0 
      E_M(NZ+2,RV,AV)= 1.D0 
      D_M(1,RV,AV)=0.D0
      D_M(NZ+2,RV,AV)=0.D0
   ENDDO
ENDDO


DO ZV=1,NZ+2
   DO AV=1,NA
      E_M(ZV,NR+2,AV)= .2D0
      D_M(ZV,NR+2,AV)= 0.D0
   ENDDO
ENDDO

EMISSAO = 0.D0

!MIST=1 (10% H2O 10% CO2)
!MIST=2 (20% H2O 10% CO2)
!MIST=3 (CO2:QUALQUER FRA��O)
!MIST=4 (H2O:QUALQUER FRA��O)
!MIST=5 (G�S CINZA: DEFINIR CA_CINZA)
MIST=5

!COEFICIENTE DE ABSOR��O PARA G�S CINZA: � ATIVADO SOMENTE QUANDO MIST=5

!---------------------------------------------------------------
! DADOS PARA PROCEDIMENTO DE CALCULO DOS COEFICIENTES DE ABSORCAO E ESPELHAMENTO DESCRITO NO PAPER.
FV  = 1.D-5 !1.D-3 !1.D-7 !11.D0/100.D0  !1.D-5 ! FRACAO VOLUMAR DAS NANOSHELLS.
R_N = 75.D-9 !10.D-9 !40.D-9 RAIO DAS NANOSHELLSSSS

! NN = FV/(4.D0/3.D0*PI*R_N**3.D0)
!WRITE(*,*) NN

N_TE = 361

! 1: Difuso, 2:HG, 3: Rayleigh, 4: Mie  
CPART = 1
CTEC = 1

!Domb = 0
	
G = 0.D0
!DESREDP = (1.D0 - 0.D0) 
!DESREDT = (1.D0 - 0.d0)


!MiePlot para nanosphere, R = 60 mm, nt =1,45 e lambda = 532 nm
Qs =  2.3d0
Qa =  2.04d0 

!MiePlot para nanoshell, R1 = 61.5 mm, R2 = 75 nm, nt =1,45 e lambda = 632.8 nm
!Qa = 2.26d0
!Qs = 2.25d0


!MiePlot para nanoshell, R1 = 14.5 mm, R2 = 20 nm, nt =1,45 e lambda = 632.8 nm
!Qa = 9.059d0
!Qs = 1.567d0

!MiePlot para nanoshell, R1 = 55 mm, R2 = 65 nm, nt =1,45 e lambda = 820 nm
!Qs =  2.79d0
!Qa =  7.d-2

!MiePlot para nanoshell, R1 = 45 mm, R2 = 50 nm, nt =1,45 e lambda = 632.8 nm
!Qs =  .155d0
!Qa =  0.277d0

!MiePlot para nanoshell, R1 = 120 mm, R2 = 155 nm, nt =1,45 e lambda = 632.8 nm
!Qs =  3.90210484010685 
!Qa =  0.649735969497625

! Mie para nanosphere, raio = 75 mm, nt =1,45 e lambda = 550 nm
!QA = 1.846D0 
!QS = 3.1033D0 


! Mie para nanosphere, raio = 40 mm, nt =1,45 e lambda = 550 nm
!QA = 3.5423D0 
!QS = 3.1268D0 


! Mie para nanosphere, raio = 50 mm, nt =1,45 e lambda = 632.8 nm 
!QA = 0.91540D0 
!QS = 4.7026D0 

! Rayleigh para nanosphere, raio = 40 mm, nt =1,45 e lambda = 550 nm 
!QA =  4.7515
!QS =  4.8026 

! Mie para nanoshell (Dombreovsky et al, 2011), raio = 20 mm, delta = 0.725, nt =1,45 e lambda = 632.8 nm
!QA = 7.828D0 
!QS = 1.144D0 

! Rayleigh para nanoshlel (Dombreovsky et al, 2011), raio = 20 nm, delta = 0,725, nt =1,45 e lambda = 0.6328
!QA =  8.1414D0
!QS =  1.3610D0 

! SET VOLUMES WITH NANOPARTICLES 
!NANO_IN_VOLUME = 0 does not means that there is not nanoparticles in the volume, but that, 
!if there are particles in the volume, they are being accounted by the atomic mix approximation.
! However, in the cases that there are not nanoparticles, NANO_IN_VOLUME must be set equal to zero.    


DO ZV = 1,NZ+2
   DO RV=1,NR+2
      DO AV = 1,NA
         NANO_IN_VOLUME(ZV,RV,AV) =  0
      ENDDO
   ENDDO
ENDDO 


!DO ZV = 129,146 !*F_MULT+1 !NZ+2
DO ZV = 1,17 !*F_MULT+1 !NZ+2  
   DO RV=1,NR+2
      DO AV = 1,NA
         NANO_IN_VOLUME(ZV,RV,AV) =  0
      ENDDO
   ENDDO
ENDDO


RN2 = R_N
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (SUM(NANO_IN_VOLUME(:,:,:)) > 0 ) THEN
      QEXT = QA + QS
      R_N = (QEXT*RN2**2.d0)**.5D0 !R_N*QEXT
      FV =  FV*(R_N**3.D0)/(RN2**3.D0)
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
R_VC = R_N  


 DO ZV=1,NZ+2
    DO RV=1,NR+2
       DO AV=1,NA
          CA_CINZA(ZV,RV,AV) =  15.64D0 !+ .75D0*FV*(QA/R_N)
          CESP(ZV,RV,AV) =  0.D0 !+ .75D0*FV*(QS/R_N)
          CEXT(ZV,RV,AV) = CA_CINZA(ZV,RV,AV) + CESP(ZV,RV,AV)
       ENDDO
   ENDDO
ENDDO


! PARAMETROS DO LASER
FWHM = 5.D-3 
P_LASER = 0.35D0 
R0=FWHM/(2.D0*LOG(2.D0))**.5D0
E0 = 2.D4 ![W/(M^2)] 


! TEMPERATURAS DAS PAREDES
DO RV=1,NR+2
   DO AV=1,NA
      T_M(1,RV,AV)= 0.D0
      T_M(NZ+2,RV,AV)= 0.D0
   ENDDO
ENDDO

T_M(1,2,1)= 1000.D0

DO ZV=1,NZ+2
   DO AV=1,NA
      T_M(ZV,NR+2,AV)=0.D0 !0.D0
      T_M(ZV,1,AV)= 0.d0
   ENDDO
ENDDO

DO RV=2,NR+1
      P_R(RV)=DIMR/(2*NR)+(DIMR/NR)*(RV-2)
      P_R(1)=0
      P_R(NR+2)=DIMR
   DO ZV=2,NZ+1
      P_Z(ZV)=DIMZ/(2*NZ)+(DIMZ/NZ)*(ZV-2)
      P_Z(NZ+2)=DIMZ
   ENDDO
ENDDO

DO AV=1,NA
   P_A(AV)=DIM_ALPHA/(2*NA) + (AV-1)*DIM_ALPHA/NA
ENDDO


P_Z(NZ+2)=DIMZ
P_R(1)=0
P_R(NR+2)=DIMR

DO RV=2,NR+1
   DO ZV=2,NZ+1
      DO AV=1,NA
         T_M(ZV,RV,AV) = 0.D0
      ENDDO
   ENDDO
ENDDO



! PRESSOES PARCIAIS NO GAS
PW=.2D0
PC=.1D0 


!***************************************************************************************************
!NUMERO DE PACOTES EMITIDOS EM CADA ZONA DE SUPERFICIE E DE MEIO PARTICIPANTE (NAO INCLUI O LASER)
!***************************************************************************************************


DO RV=2,NR+1
   DO AV=1,NA
      N_Z1_M(1,RV,AV)=0
   ENDDO
ENDDO

N_Z1_M(1,2,1)=1000000

DO RV=2,NR+1
   DO AV=1,NA
      N_Z2_M(NZ+2,RV,AV)=0
   ENDDO
ENDDO

DO ZV=2,NZ+1
   DO AV=1,NA
      N_R_M(ZV,NR+2,AV)=0
   ENDDO
ENDDO

DO ZV=2,NZ+1
   DO RV=2,NR+1
      DO AV=1,NA
         N_G_M(ZV,RV,AV)=0
      ENDDO
   ENDDO
ENDDO


WRITE(raio,'(i0)')INT(RN2)
WRITE(pindex,'(i0)')CPART
WRITE(tindex,'(i0)')CTEC


DO I_DESV = 1,N_DESV


WRITE(*,*)'I_DESV',I_DESV


CALL CILINDRO_NUCLEO(DDZ,P_Z,DDR,R_M,P_R,N_G_M,N_R_M,N_Z1_M,N_Z2_M,&
&N_LASER1,Q_LASER1,Z_LASER1,R_LASER1,AL_LASER1,GA_LASER1,BE_LASER1,A_M,V_M,EMIT_M,EAB2,&
&T_M,E_M,D_M,CA_CINZA,CESP,QA,QS,DIMR_LASER1,&
&CEXT,NANO_IN_VOLUME,CESP_NANO,CABS_NANO,CEXT_NANO,R_N,FV,R_VC,FV_VC,G,N_TE,CPART,CTEC,raio)

taxa_total(I_DESV) = 0.d0 

DO ZV = 2,NZ+1
   DO RV=2,NR+1
      DO AV=1,NA

         QG(ZV-1,RV-1,AV) = QG(ZV-1,RV-1,AV) + (EAB2(ZV,RV,AV)-EMIT_M(ZV,RV,AV))/(1000.d0*V_M(ZV,RV,AV)*N_DESV)
                   
         taxa_total(I_DESV) = taxa_total(I_DESV) + (EAB2(ZV,RV,AV)-EMIT_M(ZV,RV,AV))

      ENDDO
    ENDDO
ENDDO


DO RV=2,NR+1
   DO AV=1,NA
      QWX1(RV-1,AV) = QWX1(RV-1,AV) + (EAB2(1,RV,AV)-EMIT_M(1,RV,AV))/(1000*A_M(1,RV,AV)*N_DESV)
      QWX2(RV-1,AV) = QWX2(RV-1,AV) + (EAB2(NZ+2,RV,AV)-EMIT_M(NZ+2,RV,AV))/(1000*A_M(NZ+2,RV,AV)*N_DESV)
   ENDDO
ENDDO


DO ZV = 2,NZ+1
   DO AV=1,NA
      QWX3(ZV-1,AV) = QWX3(ZV-1,AV) + (EAB2(ZV,NR+2,AV)-EMIT_M(ZV,NR+2,AV))/(1000.d0*A_M(ZV,NR+2,AV)*N_DESV)
   ENDDO
ENDDO

DO ZV=1,NZ+2
   DO AV=1,NA
      QBG(ZV,1,AV) = 0.D0
      QBG(ZV,NR+2,AV) = 0.D0
   ENDDO
ENDDO


DO RV=2,NR+1
   DO AV=1,NA
      QBG(1,RV,AV) = - (EAB2(1,RV,AV)-EMIT_M(1,RV,AV))/(1000*A_M(1,RV,AV))
      QBG(NZ+2,RV,AV) = 0.D0
   ENDDO
ENDDO


DO ZV=2,NZ+1
   DO RV=2,NR+1
      DO AV=1,NA
         QBG(ZV,RV,AV) = QBG(ZV-1,RV,AV) - &
         &((EAB2(ZV,RV,AV)-EMIT_M(ZV,RV,AV))/(1000*V_M(ZV,RV,AV)) )*(DIMZ/NZ)
      ENDDO
   ENDDO
ENDDO

!!DO ZV=10*F_MULT+1,40*F_MULT+1
!DO ZV=2,40*F_MULT+1
DO ZV = 14,16
   DO RV=29,31
      DO AV=1,NA
         SAMPLE(I_DESV) = SAMPLE(I_DESV) + (EAB2(ZV,RV,AV)-EMIT_M(ZV,RV,AV))
      ENDDO
    ENDDO
ENDDO


ENDDO


!***********************************************************************************
!RESULTADOS  
!***********************************************************************************



!***********************************************************************************
!DADOS PARA CALCULO DO DESVIO PADRAO 
!***********************************************************************************

OPEN(3101,FILE='./resultados/taxa_total'//trim(pindex)//''//trim(tindex)//'.m')

   WRITE(3101,*)'taxa_total=['
      DO I_DESV = 1,N_DESV
         WRITE(3101,*)taxa_total(I_DESV)
      ENDDO
   WRITE(3101,*)'];'

CLOSE(3101)

!***********************************************************************************

!***********************************************************************************

OPEN(31,FILE='./resultados/radiacao'//trim(pindex)//''//trim(tindex)//'.m')

   WRITE(31,*)'NZ=',NZ,';'
   WRITE(31,*)'NR=',NR,';'
   WRITE(31,*)'NA=',NA,';'
   WRITE(31,*)'N_DESV=',N_DESV,';'

  
   WRITE(31,*)'QWX1=['
      DO RV=1,NR
         DO AV=1,NA
            WRITE(31,*)QWX1(RV,AV)
         ENDDO
      ENDDO
   WRITE(31,*)'];'

   WRITE(31,*)'QWX2=['
      DO RV=1,NR
         DO AV=1,NA
            WRITE(31,*)QWX2(RV,AV)
         ENDDO
      ENDDO
   WRITE(31,*)'];'

   WRITE(31,*)'QWX3=['
      DO ZV=1,NZ
         DO AV=1,NA
            WRITE(31,*)QWX3(ZV,AV)
         ENDDO
      ENDDO
   WRITE(31,*)'];'


  WRITE(31,*)'QG=['
      DO ZV=1,NZ
         DO RV=1,NR
            DO AV=1,NA
              WRITE(31,*)QG(ZV,RV,AV)
            ENDDO
         ENDDO
      ENDDO
   WRITE(31,*)'];'


   WRITE(31,*)'TAXA1=['
      DO RV=1,NR
         DO AV=1,NA
            WRITE(31,*)QWX1(RV,AV)*A_M(1,RV+1,AV)
         ENDDO
      ENDDO
   WRITE(31,*)'];'

   WRITE(31,*)'TAXA2=['
      DO RV=1,NR
         DO AV=1,NA
            WRITE(31,*)QWX2(RV,AV)*A_M(Nz+2,RV+1,AV)
         ENDDO
      ENDDO
   WRITE(31,*)'];'

   WRITE(31,*)'TAXA3=['
      DO ZV=1,NZ
         DO AV=1,NA
            WRITE(31,*)QWX3(ZV,AV)*A_M(ZV+1,NR+2,AV)
         ENDDO
      ENDDO
   WRITE(31,*)'];'

  
 WRITE(31,*)'TAXA4=['
      DO ZV=1,NZ
         DO RV=1,NR
            DO AV=1,NA
              WRITE(31,*)QG(ZV,RV,AV)*V_M(ZV+1,RV+1,AV)
            ENDDO
         ENDDO
      ENDDO
   WRITE(31,*)'];'

   WRITE(31,*)'EMIT1=['
      DO ZV=1,NZ+2
         DO RV=1,NR+2
            DO AV=1,NA
              WRITE(31,*)EMIT_M(ZV,RV,AV)
            ENDDO
         ENDDO
      ENDDO
   WRITE(31,*)'];'


   WRITE(31,*)'EAB2=['
      DO ZV=1,NZ+2
         DO RV=1,NR+2
            DO AV=1,NA
              WRITE(31,*)EAB2(ZV,RV,AV)
            ENDDO
         ENDDO
      ENDDO
   WRITE(31,*)'];'
   
    WRITE(31,*)'A_M1=['
      DO ZV=1,NZ+2
         DO RV=1,NR+2
            DO AV=1,NA
              WRITE(31,*)A_M(ZV,RV,AV)
            ENDDO
         ENDDO
      ENDDO
   WRITE(31,*)'];'


   WRITE(31,*)'V_M1=['
      DO ZV=1,NZ+2
         DO RV=1,NR+2
            DO AV=1,NA
              WRITE(31,*)V_M(ZV,RV,AV)
            ENDDO
         ENDDO
      ENDDO
   WRITE(31,*)'];'

   !TROCA DE CALOR NO G�S
    WRITE(31,*)'QB=['
      DO ZV=1,NZ
         DO RV=1,NR
            DO AV=1,NA
              WRITE(31,*)QBG(ZV,RV,AV)
            ENDDO
         ENDDO
      ENDDO
   WRITE(31,*)'];'

   WRITE(31,*)'NT=',0,';'

   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for RV=',1,':NR'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'qfw1(RV,AV)=QWX1(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
 
   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for RV=',1,':NR'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'qfw2(RV,AV)=QWX2(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'

   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'qfw3(ZV,AV)=QWX3(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'

   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ'
   WRITE(31,*)'for RV=',1,':NR'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'dq(ZV,RV,AV)=QG(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   WRITE(31,*)'end'


   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for RV=',1,':NR'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'qw1(RV,AV)=TAXA1(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
 
   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for RV=',1,':NR'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'qw2(RV,AV)=TAXA2(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'

   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'qw3(ZV,AV)=TAXA3(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'


   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ'
   WRITE(31,*)'for RV=',1,':NR'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'dq_taxa(ZV,RV,AV)=TAXA4(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   WRITE(31,*)'end'

   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ+2'
   WRITE(31,*)'for RV=',1,':NR+2'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'emit(ZV,RV,AV)=EMIT1(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   
   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ+2'
   WRITE(31,*)'for RV=',1,':NR+2'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'eab(ZV,RV,AV)=EAB2(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   
   
   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ+2'
   WRITE(31,*)'for RV=',1,':NR+2'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'a_m(ZV,RV,AV)=A_M1(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   
   
   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ+2'
   WRITE(31,*)'for RV=',1,':NR+2'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'v_m(ZV,RV,AV)=V_M1(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   
   WRITE(31,*)'NT=',0,';'
   WRITE(31,*)'for ZV=',1,':NZ'
   WRITE(31,*)'for RV=',1,':NR'
   WRITE(31,*)'for AV=',1,':NA'
   WRITE(31,*)'NT=NT+',1,';'
   WRITE(31,*)'q(ZV,RV,AV)=QB(NT);'
   WRITE(31,*)'end'
   WRITE(31,*)'end'
   WRITE(31,*)'end'

    WRITE(31,*)'z=['
      DO ZV=2,NZ+1
         WRITE(31,*)P_Z(ZV)
      ENDDO
   WRITE(31,*)'];'


   WRITE(31,*)'r=['
      DO RV=2,NR+1
         WRITE(31,*) P_R(RV)
      ENDDO
   WRITE(31,*)'];'


    WRITE(31,*)'alpha=['
      DO AV=2,NA
         WRITE(31,*)P_A(AV)
      ENDDO
    WRITE(31,*)'];'  

CLOSE(31)

OPEN(32,FILE='./resultados/dq'//trim(pindex)//''//trim(tindex)//'.txt')
      DO ZV=1,NZ
         DO RV=1,NR
            DO AV=1,NA
               WRITE(32,*)QG(ZV,RV,AV)*1000.D0
            ENDDO
         ENDDO
      ENDDO
CLOSE(32)

MEAN = SUM(SAMPLE)/N_DESV

SDESV = SUM((SAMPLE - MEAN)**2.D0)

DESV = (SDESV/(N_DESV-1))**.5D0

WRITE(*,*)MEAN,DESV*100.D0/MEAN

CONFIDENCE = 0.D0 

DO I_DESV =1,N_DESV

   IF (SAMPLE(I_DESV) .GE. MEAN - MEAN*0.01D0 .AND. SAMPLE(I_DESV) .LE. MEAN + MEAN*0.01d0) THEN   

      CONFIDENCE = CONFIDENCE + 1.D0

   ENDIF

ENDDO

CONFIDENCE = CONFIDENCE/N_DESV

WRITE(*,*)CONFIDENCE

OPEN(35,FILE='./resultados/standart_deviation'//trim(pindex)//''//trim(tindex)//'.m')

WRITE(35,*)'sig=',DESV*100.D0/MEAN ! percentual

WRITE(35,*)'ci=',confidence

close(35)

END PROGRAM CILINDRO_PPT

