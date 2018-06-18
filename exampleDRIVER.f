      module dataArray

c     This is a module. It holds variables that can be made available
c     anywhere. These are frequently used and I do this so I don't 
c     have to pass them all of the time. 
      DOUBLE PRECISION,allocatable :: T(:,:)
      INTEGER epochPts,numEpochs
      DOUBLE PRECISION dt

      end module dataArray

      PROGRAM driver
c     This is the test driver for the modified hierarchical 
c     decomposition method. It does the following:
c     1. Loads the data set - figure 6 in the paper
c     2. Assigns parameters 
c     3. Calls the method
c     4. Outputs the results (as a file)

      use dataArray
      use ARmatrices

      IMPLICIT NONE

      DOUBLE PRECISION,allocatable :: PC(:,:),L(:,:),HD(:,:),Lout(:,:),
     $     Tsave(:,:),weights(:,:)
      DOUBLE PRECISION,allocatable :: pvar(:)
      DOUBLE PRECISION fLo,fHi,fStep,epochLength
      INTEGER M,N,nComps,i,sampRate,fLoop
      CHARACTER*25 filename,dataOut,spaceOut
c --------------------------------------------------------
c     Required user supplied parameters
c     
      sampRate = 1 ! Sampling rate 
      numEpochs = 1 ! One epoch of 512 time points
      epochLength = 512.d0 ! This is in time units, whatever you prefer
      M = 16 ! 16 channels
      N = 512 ! Total points (1 epoch of 512 points)
      epochPts = 512 ! Points per epoch

      nComps = 3 ! Number of PCS to use
      numLags = 2 ! Number of lags to use.

      filename = 'Xcyc16.dat' ! Data file, channels are columns
      spaceOut = 'XcyclicT.dat' ! Output (spatial weights) file

      flo = 0.d0 ! Lower bound of frequency range
      fhi = 0.5d0 ! Upper bound of frequency range

      numFreqs = 11 ! Number of (evenly spaced) frequency points
      
      dt = 1.d0/DBLE(sampRate)

c     -----------------------------
c     For basic usage, one doesn't need to change anything below 
c     this
c     -----------------------------

      ! Allocate the lagList
      ALLOCATE(lagList(numLags))
      ! Allocate the spatial weights (where the results will be stored)
      ALLOCATE(weights(nComps,M))

      ! LOAD THE FILE
      IF(ALLOCATED(T)) DEALLOCATE(T)
      ALLOCATE(T(M,N))
     
c     Load the data into the array T

      CALL loadData(filename,T,M,N)

c     Compute the Principal components

      ALLOCATE(PC(M,N))
      ALLOCATE(L(M,M))
      ALLOCATE(Lout(nComps,M))
      ALLOCATE(pvar(M))
      
      CALL PCA(M,N,T,PC,L,pvar)

      DEALLOCATE(T) ! Deallocate and reallocate T
      ALLOCATE(T(nComps,N))

      T = PC(1:nComps,:)

      DEALLOCATE(PC) ! Don't need it anymore

      IF(.NOT.allocated(HD)) THEN ! Checking (should do this everywhere)
         ALLOCATE(HD(nComps,N))
      END IF
      
      ! Form the lag list
      DO i=1,numLags
         lagList(i) = i
      END DO
      
      IF(allocated(freq)) DEALLOCATE(freq)
      ALLOCATE(freq(numFreqs))
      
      fStep = (fHi-fLo)/DBLE(numFreqs-1)

      DO i=1,numFreqs
         freq(i) = DBLE(i-1)*fStep+fLo
      END DO

c     Call the FEHD algorithm.

      CALL FEHD(nComps,M,N,L,HD,Lout)
      
      
      ! Store and write the results
      weights = Lout
         
      OPEN(16,file=spaceOut,STATUS='NEW')

      DO i=1,nComps
         write(16,*) weights(i,:)
      END DO

      CLOSE(16)
      
      END

      SUBROUTINE loadData(filename,T,M,N)

      IMPLICIT NONE

      CHARACTER*25 filename
      INTEGER M,N
      DOUBLE PRECISION,dimension(M,N) :: T
      INTEGER i
      OPEN(15,FILE=filename,STATUS='OLD')

      i = 1
      
      DO 

         READ(15,FMT=*,END=100) T(:,i)
         i = i + 1

         IF(i.GT.N) GOTO 100

      END DO

 100  CONTINUE
      
      CLOSE(15)

      END SUBROUTINE
      
