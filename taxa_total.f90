! First subroutine to write taxa_total array
SUBROUTINE write_taxa_total(N_DESV, taxa_total, pindex, tindex)
  INTEGER, INTENT(IN) :: N_DESV
  REAL, INTENT(IN) :: taxa_total(N_DESV)
  CHARACTER(LEN=*), INTENT(IN) :: pindex, tindex

  INTEGER :: i_desv
  OPEN(3101, FILE='./resultados/taxa_total'//TRIM(pindex)//''//TRIM(tindex)//'.m')
  WRITE(3101, *) 'taxa_total=['
  DO i_desv = 1, N_DESV
    WRITE(3101, *) taxa_total(i_desv)
  END DO
  WRITE(3101, *) '];'
  CLOSE(3101)
END SUBROUTINE write_taxa_total