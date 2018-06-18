      subroutine pca(M,N,A,PC,TFM,pvar)
c     INPUTS
c       M,N (integer) - size of the data array A
c       A (double) - the data array
c     OUTPUTS
c       PC (double) - the principal components

      IMPLICIT NONE
      INTEGER M
      INTEGER N
      DOUBLE PRECISION A(M,N)
      DOUBLE PRECISION AT(N,M)
      DOUBLE PRECISION PC(M,N)
      DOUBLE PRECISION TFM(M,M)
      DOUBLE PRECISION pvar(M)
      INTEGER i,LWORK
      INTEGER INFO
      
      DOUBLE PRECISION U(M,M),VT(M,M),S(M),A_out(M,M),W(M,M)
      DOUBLE PRECISION,allocatable :: WORK(:)
      INTEGER,allocatable :: IWORK(:)

      LWORK=-1

      IF(allocated(IWORK)) DEALLOCATE(IWORK)
      ALLOCATE(IWORK(8*M))

c     Form the covariance matrix
      AT = TRANSPOSE(A)

      A_out = MATMUL(A,AT)

c     Call the SVD routine...
c
      ALLOCATE(WORK(10))
      call DGESDD('A',M,M,A_out,M,S,U,M,VT,M,WORK,LWORK,IWORK,INFO)

      LWORK = INT(WORK(1))

      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))

      call DGESDD('A',M,M,A_out,M,S,U,M,VT,M,WORK,LWORK,IWORK,INFO)

      W(:,:) = 0.0
      do i=1,M
         IF(S(i)>1e-3) THEN
            W(i,i) = 1.0/sqrt(S(i))
         ELSE
            W(i,i) = 1000.0
         END IF
      end do
      
      pvar = 100.0d0*S/SUM(S)

      TFM = MATMUL(W,TRANSPOSE(U))
!      call DGEMM('N','T',M,M,M,1.0d0,W,M,U,M,0.0d0,TFM,M)

c     Apply the transformation to the data

      PC = MATMUL(TFM,A)

!      call DGEMM('N','N',M,N,M,1.0d0,TFM,M,A,M,0.0d0,PC,M)

      END SUBROUTINE
