SUBROUTINE desvio_padrao(I_DESV,N_DESV,taxa_total,pindex,tindex)

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

END subroutine