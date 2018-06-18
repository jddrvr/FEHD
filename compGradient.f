      SUBROUTINE compGradient(numParticles,numDim,particle,gradient,
     $     compute)

      USE ARmatrices, only : A
      
      IMPLICIT NONE
c     ---------------------------------------
c     IN/OUT

      INTEGER, INTENT(IN) :: numParticles,numDim
      DOUBLE PRECISION,INTENT(IN),DIMENSION(numDim,numParticles) ::
     $     particle
      LOGICAL, INTENT(IN),DIMENSION(numParticles) :: compute
      DOUBLE PRECISION,INTENT(OUT),DIMENSION(numDim,numParticles) ::
     $     gradient
c     ---------------------------------------
c     Counters

      INTEGER :: particleNumber,dim
c     ---------------------------------------

c     Functions

      INTEGER :: OMP_GET_THREAD_NUM

c     data holders

      DOUBLE PRECISION,DIMENSION(numDim) :: particleL,particleR
      DOUBLE PRECISION :: causalityL,causalityR

c     Parameters

      DOUBLE PRECISION :: dt,Asave(SIZE(A,1),SIZE(A,2))
      PARAMETER(dt=1e-8)
      
      Asave = A

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(particleNumber,
!$OMP& dim,particleL,particleR,causalityL,causalityR)

      A = Asave

!$OMP DO SCHEDULE(DYNAMIC,5)
      DO particleNumber=1,numParticles

         IF(.NOT.compute(particleNumber)) CYCLE
         
         DO dim = 1,numDim
         
            particleL = particle(:,particleNumber)
            particleR = particle(:,particleNumber)

            particleL(dim) = particle(dim,particleNumber)-dt
            particleR(dim) = particle(dim,particleNumber)+dt

            CALL grangerInt(numDim,particleL,causalityL)
            CALL grangerInt(numDim,particleR,causalityR)

            gradient(dim,particleNumber) = 
     $           (causalityR-causalityL)/(2.0*dt)

         END DO

      END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL

      END SUBROUTINE
