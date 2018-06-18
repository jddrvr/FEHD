      SUBROUTINE calcDet(det,M,A)

      IMPLICIT NONE


      INTEGER M
      COMPLEX*16,dimension(M,M) :: A,L,U
      DOUBLE PRECISION :: det
      INTEGER j,k

      U = A
      L = 0

      DO k=1,M-1
         
         DO j=k+1,M

            L(j,k) = U(j,k)/U(k,k)
            U(j,k:M) = U(j,k:M) - L(j,k)*U(k,k:M)

         END DO

      END DO

      det = 1.0d0

      DO k=1,M
         det = det*DBLE(U(k,k))
      END DO

      END SUBROUTINE
      


      
