      SUBROUTINE FEHD(nComps,M,N,L,HD,Lout)
      
      use ARmatrices
      use dataArray
      IMPLICIT NONE
      
      INTEGER nComps,M,N
      DOUBLE PRECISION,dimension(nComps,N) :: HD,Ttmp
      DOUBLE PRECISION,dimension(M,M) :: L
      DOUBLE PRECISION,dimension(nComps,M) :: Lout,Lstep
      DOUBLE PRECISION,dimension(M) :: pvar
      DOUBLE PRECISION,allocatable :: Rdecor(:,:)
      DOUBLE PRECISION,allocatable :: decor(:,:),TMP(:,:),
     $     invdecor(:,:),Q(:,:),Qtmp(:,:),Rtrans(:,:)
      DOUBLE PRECISION :: GI
      INTEGER Rpoints,i,j,dim,numParticles
      DOUBLE PRECISION,allocatable :: theta(:),blo(:),bup(:)

      numEpochs = N/epochPts

      numParticles = 5000 ! This sets the number of particles
      
      Lstep = L(1:nComps,:)

      Lout = 0.d0
      
      DO comps = nComps,2,-1

         dim = comps-1

         IF(allocated(theta)) DEALLOCATE(theta)
         ALLOCATE(theta(dim))

         IF(allocated(bup)) DEALLOCATE(bup)
         ALLOCATE(bup(dim))

         IF(allocated(blo)) DEALLOCATE(blo)
         ALLOCATE(blo(dim))

         blo = -3.1415/2.0
         bup = 3.1415/2.0
         
c     Construct the auto-regressive model.
         
         CALL makeAR()

c     Make the residuals ortho-normal.

         Rpoints = size(R,1)

         IF(allocated(Rdecor)) DEALLOCATE(Rdecor)
         ALLOCATE(Rdecor(comps,Rpoints))
         IF(allocated(TMP)) DEALLOCATE(TMP)
         ALLOCATE(TMP(comps,comps))
         IF(allocated(decor)) DEALLOCATE(decor)
         IF(allocated(invdecor)) DEALLOCATE(invdecor)
         IF(allocated(Q)) DEALLOCATE(Q)
         IF(allocated(Qtmp)) DEALLOCATE(Qtmp)
         IF(allocated(Rtrans)) DEALLOCATE(Rtrans)
         ALLOCATE(Qtmp(comps,comps))
         ALLOCATE(Q(comps,comps))
         ALLOCATE(decor(comps,comps))
         ALLOCATE(invdecor(comps,comps))
         ALLOCATE(Rtrans(comps,Rpoints))

         Rtrans = TRANSPOSE(R)
         
         CALL pca(comps,Rpoints,Rtrans,Rdecor,decor,pvar)

         R = Rdecor

         invdecor = decor

         CALL invert(invdecor,comps)

         DO i=1,numLags
            TMP = A((i-1)*comps+2:i*comps+1,:)
            TMP = MATMUL(TRANSPOSE(invdecor),TMP)
            TMP = MATMUL(TMP,TRANSPOSE(decor))
            A((i-1)*comps+2:i*comps+1,:)= TMP
         END DO

         ! Find the rotations using steepest descent. 
         CALL STEEPESTDESCENT(theta,blo,bup,dim,numParticles)

         ! Rotate the data
         Q = 0.0
         
         DO i=1,comps
            Q(i,i) = 1.0d0
         END DO
         
         DO i=1,dim
               
            Qtmp = 0.0

            DO j=1,comps
               Qtmp(j,j) = 1.0d0
            END DO
               
            Qtmp(i,i) = COS(theta(i))
            Qtmp(i,comps) = -SIN(theta(i))
            Qtmp(comps,i) = SIN(theta(i))
            Qtmp(comps,comps) = COS(theta(i))

            Q = MATMUL(Qtmp,Q)
         END DO
         
         TMP = MATMUL(Q,decor)

         Lstep(1:comps,:) = MATMUL(TMP,Lstep(1:comps,:))

         Lout(comps,:) = Lstep(comps,:)
         
!     Left multiply to get the data.

         Ttmp(1:comps,:) = MATMUL(TMP,T)

         HD(comps,:) = Ttmp(comps,:)

         DEALLOCATE(T)
         ALLOCATE(T(comps-1,N))

         T = Ttmp(1:comps-1,:)

         write(*,*) 'Component finished'

      END DO

      HD(1,:) = T(1,:)

      Lout(1,:) = Lstep(1,:)

      IF(allocated(Rdecor)) DEALLOCATE(Rdecor)
      IF(allocated(TMP)) DEALLOCATE(TMP)
      IF(allocated(decor)) DEALLOCATE(decor)
      IF(allocated(invdecor)) DEALLOCATE(invdecor)
      IF(allocated(Q)) DEALLOCATE(Q)
      IF(allocated(Qtmp)) DEALLOCATE(Qtmp)


      END
