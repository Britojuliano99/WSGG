SUBROUTINE resultados(pindex,tindex)


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

end SUBROUTINE