      subroutine ARspectrum()

      use ARmatrices,only:A,spec,Tf,numFreqs,numLags,comps,
     $     lagList,freq
      use dataArray,only:dt
      
      IMPLICIT NONE

      DOUBLE PRECISION pi
      INTEGER i,j
      COMPLEX*16 :: argmt

      ! Fortran doesn't know this

      pi = 4.d0*ATAN(1.d0)

      ! Allocate and initialize some things

      spec = DCMPLX(0.d0,0.d0)
         
      Tf = DCMPLX(0.d0,0.d0)
      
      DO i=1,numFreqs

         DO j=1,comps
            Tf(j,j,i) = DCMPLX(1.d0,0.d0)
         END DO

         DO j=1,numLags

            argmt = -2.d0*pi*DBLE(lagList(j))*dt*freq(i)*
     $           DCMPLX(0.d0,1.d0)

            Tf(:,:,i) = Tf(:,:,i)-
     $           TRANSPOSE(A(2+(j-1)*comps:j*comps+1,:))*
     $           EXP(argmt)
           
         END DO

         spec(1:comps,1:comps,i) =
     $        MATMUL(TRANSPOSE(CONJG(Tf(:,:,i))),Tf(:,:,i))

      END DO
      
      END SUBROUTINE
      
