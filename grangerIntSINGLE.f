      SUBROUTINE grangerInt(dim,particle,causality)

      use ARmatrices,only:A,spec,Tf,numFreqs,numLags,freq,
     $     comps
      
      IMPLICIT NONE

      INTEGER dim
      DOUBLE PRECISION :: causality
      DOUBLE PRECISION,dimension(dim) :: particle
       
      INTEGER i,j,k
      COMPLEX*16,dimension(dim,dim) :: Hxx,Sself
      COMPLEX*16,dimension(dim,dim) :: specSub
      COMPLEX*16,dimension(dim,1) :: column
      COMPLEX*16,dimension(1,dim) :: row
      DOUBLE PRECISION gIndex(numFreqs)
      DOUBLE PRECISION,dimension(2,comps) :: Arows
      DOUBLE PRECISION,dimension(comps*numLags+1,comps) :: Asave
      DOUBLE PRECISION,dimension(comps,comps) :: Atmp
      DOUBLE PRECISION detSelf,detFull,cc,ss
      DOUBLE PRECISION,dimension(comps,comps) :: Q,Qtmp

c     Rotate A here, and pass to ARspectrum
      
      Asave = A

      IF(ALLOCATED(spec)) DEALLOCATE(spec)
      ALLOCATE(spec(comps,comps,numFreqs))
      IF(ALLOCATED(Tf)) DEALLOCATE(Tf)
      ALLOCATE(Tf(comps,comps,numFreqs))

      A = Asave

      Q = 0.0

      DO j=1,comps
         Q(j,j) = 1.0d0
      END DO
         
      DO j=1,dim

         Qtmp = 0.0
         DO k=1,comps
            Qtmp(k,k) = 1.0d0
         END DO
         
         cc = COS(particle(j))
         ss = SIN(particle(j))
         
         Qtmp(j,j) = cc
         Qtmp(j,comps) = -ss
         Qtmp(comps,j) = ss
         Qtmp(comps,comps) = cc
         
         Q = MATMUL(Qtmp,Q)
         
      END DO
      
      DO j=1,numLags
         
         Atmp = A((j-1)*comps+2:j*comps+1,:)
         Atmp = MATMUL(Atmp,TRANSPOSE(Q))
         Atmp = MATMUL(Q,Atmp)
         A((j-1)*comps+2:j*comps+1,:) = Atmp
      END DO

      CALL ARspectrum()
         
      Sself = (0.d0,0.d0)

      DO j=1,numFreqs

         DO k=1,dim
            column(k,1) = Tf(k,comps,j)
            row(1,k) = Tf(comps,k,j)
         END DO

         Hxx = Tf(1:comps-1,1:comps-1,j)-
     $        (1.0,0.0)/Tf(comps,comps,j)*
     $        MATMUL(column,row)
         
         
         Sself = MATMUL(TRANSPOSE(CONJG(Hxx)),Hxx)
         
         DO k=1,dim
            column(k,1) = spec(k,comps,j)
            row(1,k) = spec(comps,k,j)
         END DO
         
         specSub = spec(1:dim,1:dim,j)-
     $        (1.0,0.0)/spec(comps,comps,j)*
     $        MATMUL(column,row)
         
         
         CALL calcDet(detSelf,comps-1,Sself)
         CALL calcDet(detFull,comps-1,specSub)
         
         gIndex(j)=LOG(detSelf/detFull)
         
      END DO

c     Perhaps a more sophisticated Quad method
      
      causality = SUM(gIndex)
      
c     Restore A for the next trial.
      
      A = Asave
      
      END SUBROUTINE
