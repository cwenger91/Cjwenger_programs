PROGRAM ODE_Solver
     IMPLICIT NONE
     REAL(KIND=8), ALLOCATABLE :: YSOL(:,:), XSOL(:), k1(:), k2(:), k3(:), k4(:), Y(:), y0(:), f(:), ym(:)
     REAL(KIND=8) :: delx, x0,x,temp
     INTEGER :: i,j,n=0,h,ierror = 0 !ierror set to 0 for success  
   !User inputs go here
    WRITE(*,*) "Step Size:"
    READ(*,*) delx
    WRITE(*,*) "Number of steps:"
    READ(*,*) h

    OPEN(unit=1,file="input.dat",Status='OLD',Action='READ',Iostat=ierror)
         !gets the size of the matrix from the file
         DO
         READ(1,*,Iostat=ierror)temp
         IF(ierror.NE.0)EXIT !Exits loop if there is an error in file	     
         n = n + 1
         END DO

    !allocates the meomery based on the size of the file
    ALLOCATE(YSOL(n,h),XSOL(h),k1(n),k2(n),k3(n),k4(n),Y(n),y0(n),f(n),ym(n))
        REWIND(Unit=1) !resets the position of reading the file
       !read intial x value
        READ(1,*) XSOL(1)
        DO i=1,n-1
        !read intial y values
            READ(1,*) YSOL(i,1)
        END DO
    CLOSE(1)
    


    x = XSOL(1)
    DO i = 1, h
        CALL RK4(n,delx,x,y,k1,k2,k3,k4)
        XSOL(i) = x !stores the value of x for the given equations
        do j=1, n
            YSOL(j,i) = Y(j) !stores the solution of the equations given
        END DO
    END DO
    !Write the solutions to a dat file, Matlab can read and plot later
    OPEN(unit=2,file="output.dat")
 !       j  YSOL(1)   YSOL(2)   YSOL(3)   YSOL(4)   XSOL
    DO j = 1, h
        WRITE(2,*)j,YSOL(1,j), YSOL(2,j), YSOL(3,j), YSOL(4,j),XSOL(j)
    END DO

    CLOSE(1)
    CLOSE(2)

END PROGRAM
!**************RK4 Subroutine*******************!
SUBROUTINE RK4(n,delx,x,y,k1,k2,k3,k4)
    IMPLICIT NONE
    REAL*8 :: Y(n),k1(n),k2(n),k3(n),k4(n),xm,delx,x,ym(n)
    INTEGER :: i,n
    
CALL Derivatives(x,y,k1,n)
    xm = x + .50d0*delx
    DO i=1,n
        ym(i)=y(i)+k1(i)*delx*.50d0
    END DO
    
CALL Derivatives(xm,ym,k2,n)
    DO i=1,n
        ym(i) = y(i) + delx*.50d0*k2(i)
    END DO
CALL Derivatives(xm,ym,k3,n)
    xm = x + delx
    DO i=1,n
        ym(i) = y(i)+k3(i)*delx
    END DO
CALL Derivatives(xm,ym,k4,n)
    DO i=1,n
        y(i) = y(i)+((delx/6.0d0)*(k1(i)+2.0d0*k2(i)+2.0d0*k3(i)+k4(i)))
    END DO
    x = xm
END SUBROUTINE
!*********************Derivatives Subroutine*******!
SUBROUTINE Derivatives(x,y,f,n)
    IMPLICIT NONE
    INTEGER :: n
    REAL*8 :: x,y(n),f(n)

    !Set order for equations
    !y(1) = y1
    !y(2) = dy1/dx
    !y(3) = y2
    !y(4) = dy2/dx
    !f(n)...
    !Derivatives listed below  
    !Test case 
    f(1) = y(2)
    f(2) = y(3)
    f(3) = SIN(x)-ABS(y(2)*y(2))*y(3)-6*y(1)*y(1)*y(1)*y(2)-y(1)*SIN(y(1))
END SUBROUTINE
