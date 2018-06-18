      module ARmatrices

      DOUBLE PRECISION,allocatable :: A(:,:)
      DOUBLE PRECISION,allocatable :: R(:,:)
      COMPLEX(8),allocatable :: spec(:,:,:)
      COMPLEX(8),allocatable :: Tf(:,:,:)
      INTEGER,allocatable :: lagList(:)
      DOUBLE PRECISION,allocatable :: freq(:)
      INTEGER numFreqs,numLags,comps

!$OMP THREADPRIVATE(A)
!$OMP THREADPRIVATE(R)
!$OMP THREADPRIVATE(spec)
!$OMP THREADPRIVATE(Tf)

      END module

      SUBROUTINE makeAR()

      use dataArray
      use ARmatrices

      IMPLICIT NONE
      
      INTEGER :: compEpochPts,usedTimes
      INTEGER :: j,k
      INTEGER :: maxLag,numCh,covSize                         
      DOUBLE PRECISION,allocatable :: LHSData(:,:),lagData(:,:),
     $     lagCov(:,:),lagLHDcov(:,:)

      
      numEpochs = size(T,2)/epochPts
      numCh = size(T,1)

      maxLag = MAXVAL(lagList,1)

      compEpochPts = epochPts - maxLag
      usedTimes = compEpochPts * numEpochs

      covSize = numCh*numLags+1
      
      IF(allocated(LHSdata)) DEALLOCATE(LHSdata)
      IF(allocated(lagData)) DEALLOCATE(lagData)
      IF(allocated(A)) DEALLOCATE(A)
      IF(allocated(lagCov)) DEALLOCATE(lagCov)
      IF(allocated(lagLHDcov)) DEALLOCATE(lagLHDcov)

      ALLOCATE(LHSdata(usedTimes,numCh))
      ALLOCATE(lagData(usedTimes,numCh*numLags+1))

      lagData(:,1) = 1.0

      DO j=1,numEpochs

         LHSdata((j-1)*compEpochPts+1:j*compEpochPts,:)=
     $        TRANSPOSE(T(:,(j-1)*epochPts+maxLag+1:j*epochPts))

         DO k=1,numLags
            
            lagData((j-1)*compEpochPts+1:j*compEpochPts,(k-1)*
     $           numCh+2
     $           :k*numCh+1) = 
     $           TRANSPOSE(T(:,(j-1)*epochPts+maxLag+1-
     $           lagList(k):
     $           j*epochPts-lagList(k)))
         END DO
         
      END DO

      ALLOCATE(A(covSize,numCh))
      
      ALLOCATE(lagCov(covSize,covSize))
      ALLOCATE(lagLHDcov(covSize,numCh))
            
c     COMPUTE THE LEAST SQUARES SOLUTION
      
      lagCov = MATMUL(TRANSPOSE(lagData),lagData)

      CALL invert(lagCov,covSize)

      lagLHDcov = MATMUL(TRANSPOSE(lagData),LHSdata)
      
      A = MATMUL(lagCov,lagLHDcov)
      
c     COMPUTE THE RESIDUAL
      
      LHSdata = LHSdata - MATMUL(lagData,A)

      IF(allocated(R)) DEALLOCATE(R)
      ALLOCATE(R(usedTimes,numCh))

      R = LHSdata

      DEALLOCATE(LHSdata)
      DEALLOCATE(lagData)
      DEALLOCATE(lagCov)
      DEALLOCATE(lagLHDcov)
      
      END SUBROUTINE
