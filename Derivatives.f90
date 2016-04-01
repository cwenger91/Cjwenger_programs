SUBROUTINE Derivatives(x,y,f,n)
    IMPLICIT NONE
    INTEGER :: n
    REAL*8 :: x,y(n),f(n)
    !y(1) = y
    !y(2) = dy/dx
    !y(3) = z
    !y(4) = dz/dx
	
    !Derivatives listed below  
    f(1) = y(2)
    f(2) = -6*y(1) + 4*y(3)
    f(3) = y(4)
    f(4) = -12*y(3) + 4*y(1)
	
	!Test case 
    !f(1) = y(2)
	!f(2) = 5.0d0*sin(2.0d0*x) - 4*abs(y(4))*y(2) - 8*(y(3)*y(3)*y(3))*y(1)
    !f(3) = y(4)
    !f(4) = -8.0d0*abs(y(2))*y(4) - 8.0d0*(y(1)*y(1)*y(1))*y(3)
    !f(n)=...

END SUBROUTINE
