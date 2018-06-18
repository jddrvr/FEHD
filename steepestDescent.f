      SUBROUTINE STEEPESTDESCENT(theta,blo,bup,dim,numParticles)

      use ARmatrices, only: A

      IMPLICIT NONE

      INTEGER dim,numParticles
      DOUBLE PRECISION,dimension(dim) :: theta,blo,bup
      DOUBLE PRECISION,dimension(dim,numParticles) :: particle,grad,
     $     particleOLD,gradOLD,particleDiff,
     $     gradDiff
      DOUBLE PRECISION,dimension(numParticles) :: gValues,gValuesOLD
      DOUBLE PRECISION,dimension(numParticles) :: gamma,norm
      DOUBLE PRECISION :: dt,TOL,minValOLD,compValue
      INTEGER partNum,angle,iteration
      INTEGER,dimension(1) :: index
      LOGICAL,dimension(numParticles) :: compute
      INTEGER,dimension(numParticles) :: sumArray
c     Initialize the particles to random locations throughout the space
      DOUBLE PRECISION,dimension(size(A,1),size(A,2)) :: Asave
      minValOLD = 10000.0
      
      dt = 1e-8
      TOL = 1e-8
      particle(1,1) = RAND(2424) ! Just getting RAND started.

      compute = .TRUE.

      DO partNum=1,numParticles
         DO angle=1,dim
            particle(angle,partNum) = (bup(angle)-blo(angle))*RAND()+
     $           blo(angle)
         END DO

      END DO

      Asave = A

c     I USE OPENMP TO PARALLELIZE A LOT OF THESE CALCULATIONS.
c     I DONT CLAIM THAT THIS IS THE BEST WAY TO DO IT. 

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(partNum)

      A = Asave

c     Go through all of the particles to get the value at iteration 0.

!$OMP DO SCHEDULE(STATIC,5)
      DO partNum=1,numParticles

         CALL grangerInt(dim,particle(:,partNum),gValues(partNum))

      END DO
!$OMP ENDDO
!$OMP END PARALLEL

      grad = 0.0
      
      DO iteration = 1,50      ! 50 is an arbitrary number

c     Compute the gradient at each particle

         gradOLD = grad

         gValuesOLD = gValues

  
c     COMPUTE THE GRADIENT OF EACH 'PARTICLE'
         CALL compGradient(numParticles,dim,particle,grad,compute)


c     THIS PARALLEL SECTION UPDATES EACH PARTICLE.
c     THERE IS A TRAP FOR THOSE THAT "OVERSHOOT", 
c     ENSURING THAT THE OBJECTIVE VALUE DECREASES ON EACH
c     ITERATION.
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(partNum)
!$OMP DO SCHEDULE(DYNAMIC,100)
         DO partNum = 1,numParticles            

            IF(.NOT.compute(partNum)) CYCLE
            
            IF(iteration.LE.1) THEN ! Getting started.
               gamma(partNum) = 1e-3

            ELSE

               gamma(partNum) = DOT_PRODUCT(particle(:,partNum)-
     $              particleOLD(:,partNum),
     $              grad(:,partNum)-gradOLD(:,partNum))


               norm(partNum) = DOT_PRODUCT(grad(:,partNum)-
     $              gradOLD(:,partNum),
     $              grad(:,partNum)-gradOLD(:,partNum))

               IF(norm(partNum).LT.TOL) THEN
                  norm(partNum) = TOL
               END IF

               IF(gamma(partNum).LT.0.0) THEN
                  gamma(partNum) = 1e-4
               ELSE
                  gamma(partNum) = gamma(partNum)/norm(partNum)
               END IF
               
            END IF
             
            particleOLD(:,partNum) = particle(:,partNum)
            
            particle(:,partNum)=particle(:,partNum)-
     $           gamma(partNum)*grad(:,partNum)

            CALL grangerInt(dim,particle(:,partNum),gValues(partNum))

            ! In the event of an "overshoot"
            DO WHILE(gValues(partNum).GT.gValuesOLD(partNum))

               gamma(partNum) = 0.5*gamma(partNum)

               particle(:,partNum)=particleOLD(:,partNum)-
     $              gamma(partNum)*grad(:,partNum)

               CALL grangerInt(dim,particle(:,partNum),gValues(partNum))

            END DO
               
         END DO
!$OMP ENDDO NOWAIT
!$OMP END PARALLEL
         
         write(*,*) iteration,MINLOC(gValues),MINVAL(gValues)

         sumArray = 0

         WHERE(compute) sumArray=1

         IF(SUM(sumArray).EQ.0) EXIT
         
      END DO

      ! The minimum particle number and angle list. 
      index = MINLOC(gValues)
      theta = particle(:,index(1))
            
      END SUBROUTINE



            
