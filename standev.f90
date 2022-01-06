!     Kuntal Ghosh
!     October 2020      

!-------------------------------------------------------------------------

      PROGRAM standard_deviation

         IMPLICIT NONE

         INTEGER*8 :: ncount,i,ios
         REAL*8    :: t,xk2,vk2,xk4,vk4,xk2dev,vk2dev,&
                      xk4dev,vk4dev,sumxk2,sumxk4,sumvk2,&
                      sumvk4,rmsxk2,rmsvk2,rmsxk4,rmsvk4,xa,va

         OPEN (1, FILE='analytical.dat',STATUS='OLD')
         OPEN (2, FILE='rk2_spring.dat',STATUS='OLD')
         OPEN (3, FILE='rk4_spring.dat',STATUS='OLD')

         ncount = 0
         no_of_lines: DO
            READ (1,*,iostat=ios)
            IF (ios/=0) EXIT no_of_lines
            ncount = ncount + 1
         END DO no_of_lines

         WRITE (*,*) "No of data points =", ncount

         REWIND(1)

         sumxk2 = 0.0d0
         sumvk2 = 0.0d0
         sumxk4 = 0.0d0
         sumvk4 = 0.0d0

         DO i = 1,ncount
            READ (1,*) t,xa,va
            READ (2,*) t,xk2,vk2
            READ (3,*) t,xk4,vk4
            xk2dev = (xk2-xa)**2
            xk4dev = (xk4-xa)**2
            vk2dev = (vk2-va)**2
            vk4dev = (vk4-va)**2           
            sumxk2 = sumxk2 + xk2dev 
            sumvk2 = sumvk2 + vk2dev
            sumxk4 = sumxk4 + xk4dev
            sumvk4 = sumvk4 + vk4dev
         END DO

         rmsxk2 = DSQRT(sumxk2/DFLOAT(ncount))
         rmsxk4 = DSQRT(sumxk4/DFLOAT(ncount))
         rmsvk2 = DSQRT(sumvk2/DFLOAT(ncount))
         rmsvk4 = DSQRT(sumvk4/DFLOAT(ncount))

         WRITE (*,*) "RMS deviation values from analytical function & 
                      are as follows ="
         WRITE (*,"(A,F16.6)") "For position (RK2) =",rmsxk2 
         WRITE (*,"(A,F16.6)") "For position (RK4) =",rmsxk4 
         WRITE (*,"(A,F16.6)") "For velocity (RK2) =",rmsvk2 
         WRITE (*,"(A,F16.6)") "For velocity (RK4) =",rmsvk4 

      END PROGRAM standard_deviation
