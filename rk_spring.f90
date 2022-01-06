      PROGRAM runge_kutta

!-------------------------------------------------------------

!     Kuntal Ghosh
!     October 2020

!-------------------------------------------------------------

         IMPLICIT NONE

         INTEGER(8) :: i,n,order
         REAL(8) :: t,x,f,v,dxdt,dvdt,pe,ke,te, &
                    k1v,k2v,k3v,k4v, &
                    k1r,k2r,k3r,k4r, &
                    k1,k2,k3,k4
         REAL(8), PARAMETER :: t0 = 0.0d0,  &
                               x0 = 0.0d0,  &
                               v0 = -1.0d0, &
                               tmax=500.0d0,&
                               dt = 0.1d0,  &
                               k = 1.0d0,   &
                               m = 1.0d0

         OPEN (2, FILE='analytical.dat',STATUS='UNKNOWN')

         PRINT *, "Enter the order of the Runge Kutta scheme ="
         READ (*,*) order

         n = INT((tmax-t0)/dt)

         t = t0
         x = x0
         v = v0
         CALL analytical

         t = t0
         x = x0
         v = v0
         rk_order: SELECT CASE (order)
         CASE(2)
            CALL rk2
         CASE(4)
            CALL rk4
         END SELECT rk_order

         CONTAINS

         SUBROUTINE rk2 

            IMPLICIT NONE

            DO i = 1,n
               ke = 0.5d0*m*v**2
               pe = 0.5d0*k*x**2
               te = ke + pe
               OPEN (1, FILE='rk2_spring.dat',STATUS='UNKNOWN')
               WRITE (1,"(6F16.6)") t,x,v,ke,pe,te

               k1v = dt*slope(1,x)
               k1r = dt*slope(0,v)
               k2v = dt*slope(1,x+k1r/2.0d0)
               k2r = dt*slope(0,v+k1v/2.0d0)

               v = v + k2v
               x = x + k2r

               t = t + dt
            END DO
        
         END SUBROUTINE rk2

         SUBROUTINE rk4 

            IMPLICIT NONE

            DO i = 1,n
               ke = 0.5d0*m*v**2
               pe = 0.5d0*k*x**2
               te = ke + pe
               OPEN (1, FILE='rk4_spring.dat',STATUS='UNKNOWN')
               WRITE (1,"(6F16.6)") t,x,v,ke,pe,te

               k1v = dt*slope(1,x)
               k1r = dt*slope(0,v)
               k2v = dt*slope(1,x+k1r/2.0d0)
               k2r = dt*slope(0,v+k1v/2.0d0)
               k3v = dt*slope(1,x+k2r/2.0d0)
               k3r = dt*slope(0,v+k2v/2.0d0)
               k4v = dt*slope(1,x+k3r)
               k4r = dt*slope(0,v+k3v)

               v = v + (k1v+k4v)/6.0d0 + (k2v+k3v)/3.0d0          
               x = x + (k1r+k4r)/6.0d0 + (k2r+k3r)/3.0d0

               t = t + dt
            END DO
        
         END SUBROUTINE rk4

         REAL(8) FUNCTION slope(aa,xv) RESULT(deriv)

            IMPLICIT NONE
            REAL(8),INTENT(IN) :: xv
            INTEGER,INTENT(IN) :: aa

            IF (aa==1) THEN
               deriv = (-k/m)*xv
            ELSE
               deriv = xv
            END IF

         END FUNCTION slope

         SUBROUTINE analytical

            IMPLICIT NONE

            REAL(8) :: t,x,v,w

            w = DSQRT(k/m)

            DO i = 1,n
               x = (-1.0d0/w)*DSIN(w*t)
               v = -DCOS(w*t)
               t = t + dt
               WRITE (2,"(3F16.6)") t,x,v
            END DO

         END SUBROUTINE analytical

      END PROGRAM runge_kutta
