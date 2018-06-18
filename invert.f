      SUBROUTINE invert(A,M)
      
      INTEGER M
      DOUBLE PRECISION A(M,M)
      INTEGER IPIV(M)
      INTEGER LWORK,INFO
      DOUBLE PRECISION WORK(2*M)

      LWORK=2*M

      CALL DGETRF(M,M,A,M,IPIV,INFO)
      IF(INFO.NE.0) THEN
         write(*,*) 'DGETRF INFO=',INFO
         write(*,*) 'Stopping'
         STOP
      END IF

      CALL DGETRI(M,A,M,IPIV,WORK,LWORK,INFO)
      IF(INFO.NE.0) THEN
         write(*,*) 'DGETRI INFO=',INFO
         write(*,*) 'Stopping'
         STOP
      END IF

      END SUBROUTINE
