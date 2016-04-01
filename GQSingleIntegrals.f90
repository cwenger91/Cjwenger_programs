      !Christopher Wenger
      !Aere 361 Fall 2015
      !Gauss Quandrature single intergral method
      !This program uses the Gauss Quandrature therom to preform numerical intergrations

      !***********************************MAIN SECTION*****************************************!
      PROGRAM GQSingle
          IMPLICIT NONE
          REAL(KIND=4) :: a,b,degreeofaccuracy,eps,integral(2) ! b-> upper bound, a -> lower bound of intergral
               
      
          WRITE(*,*) "Enter in degree of accuracy i.e X.X"
          READ(*,*)degreeofaccuracy
          WRITE(*,*) "Enter in upper bound i.e X.X"
          READ(*,*)b
          WRITE(*,*) "Enter in lower bound i.e X.X"
          READ(*,*)a
		  eps = 1.0/(10.0**degreeofaccuracy)
          CALL integrator(a,b,eps,integral)
      
          WRITE(*,*) "The integral value = ", integral(2)    
      END PROGRAM GQSINGLE !ENDING OF MAIN SECTION
!*********************************INTERGRATIOR SUBROUTINE**********************************!
     !This subroutine handles all of the integration processing
      
      SUBROUTINE integrator(a,b,eps,integral)
          IMPLICIT NONE
          REAL(KIND=4) :: x(2),b,a,c,m,eps,integral(2),fun_val(2)
          REAL(KIND=4), DIMENSION(10,10) :: t,w
          INTEGER :: counter,converge, n ! n -> # number of GAUSS points
          CALL tw_vals(t,w) 
        
          c = 0.5 * (b+a)
          m = 0.5 * (b-a)

          DO n = 2,10 
             DO counter=1,n    
             x(1) = c+m*t(n-1,counter-1)
             CALL fun(fun_val(1),x(1))
             integral(1) = integral(1) + m*w(n-1,counter-1)*fun_val(1)

             x(2) = c+m*t(n,counter)
             CALL fun(fun_val(2),x(2))
             integral(2) = integral(2) + m*w(n,counter)*fun_val(2)
       
            END DO
          CALL errorCheck(eps,converge,integral)
          IF(converge .EQ. 1) THEN
          EXIT
          ELSE
            integral(1) = 0.0 !reset
            integral(2) = 0.0 !reset
          END IF
          END DO    

      END SUBROUTINE integrator !ENDING OF INTEGRATIOR SUBROUTINE
      
     !****************************************ERROR CHECKER SUBROUTINE****************************!
     !This subroutine checks to see if the value of the integration is within the range of the 
     !user's preference

      SUBROUTINE errorCheck(eps,converge,I)
        REAL(KIND=4) :: eps,I(2), Rel_error
        INTEGER :: converge
        converge = 0
        IF (I(2) .NE. 0) THEN
           Rel_error = ABS((I(2)-I(1))/I(2))
        ELSE 
           Rel_error = ABS((I(2)-I(1)/I(1)))
        END IF
    
        IF ( Rel_error .LE. eps)  THEN
           converge = 1
        END IF

      END SUBROUTINE errorCheck !ENDING OF ERROR CHCKER SUBROUTINE

     !*****************************************FUNCTION SUBROUTINE********************************!
     !This is where the user can define a new function
      SUBROUTINE fun(f,x)
          REAL(KIND=4) :: f,x
          f = (x*x)*log(2+3*x)
      END SUBROUTINE fun !ENDING OF FUCTION SUBROUTINE
            
     !****************************************TW_vals SUBROUTINE**********************************!
     !This SUBROUTINE contains the Gauss values
     !@AURTHOR Dr. Ambar K. Mitra
     !@LINK http://www.public.iastate.edu/~akmitra/aero361/design_web/g_weights.txt
      SUBROUTINE tw_vals(t,w)
      IMPLICIT NONE
      REAL(KIND=4),DIMENSION(10,10) :: t,w
      
        t=0
        t(2,1)=0.57735027
        t(2,2)=-0.57735027
        t(3,2)=0.77459667
        t(3,3)=-0.77459667
        t(4,1)=0.33998104
        t(4,2)=-0.33998104
        t(4,3)=0.86113631
        t(4,4)=-0.86113631
        t(5,2)=0.53846931
        t(5,3)=-0.53846931
        t(5,4)=0.90617985
        t(5,5)=-0.90617985
        t(6,1)=0.23861918
        t(6,2)=-0.23861918
        t(6,3)=0.66120939
        t(6,4)=-0.66120939
        t(6,5)=0.93246951
        t(6,6)=-0.93246951
        t(7,2)=0.40584515
        t(7,3)=-0.40584515
        t(7,4)=0.74153119
        t(7,5)=-0.74153119
        t(7,6)=0.94910791
        t(7,7)=-0.94910791
        t(8,1)=0.18343464
        t(8,2)=-0.18343464
        t(8,3)=0.52553241
        t(8,4)=-0.52553241
        t(8,5)=0.79666648
        t(8,6)=-0.79666648
        t(8,7)=0.96028986
        t(8,8)=-0.96028986
        t(9,2)=0.32425342
        t(9,3)=-0.32425342
        t(9,4)=0.61337143
        t(9,5)=-0.61337143
        t(9,6)=0.83603111
        t(9,7)=-0.83603111
        t(9,8)=0.96816024
        t(9,9)=-0.96816024
        t(10,1)=0.14887434
        t(10,2)=-0.14887434
        t(10,3)=0.43339539
        t(10,4)=-0.43339539
        t(10,5)=0.67940957
        t(10,6)=-0.67940957
        t(10,7)=0.86506337
        t(10,8)=-0.86506337
        t(10,9)=0.97390653
        t(10,10)=-0.97390653


        w=0
        w(2,1)=1
        w(2,2)=1
        w(3,1)=0.88888889
        w(3,2)=0.55555555
        w(3,3)=0.55555555
        w(4,1)=0.65214515
        w(4,2)=0.65214515
        w(4,3)=0.34785485
        w(4,4)=0.34785485
        w(5,1)=0.56888889
        w(5,2)=0.47862867
        w(5,3)=0.47862867
        w(5,4)=0.23692689
        w(5,5)=0.23692689
        w(6,1)=0.46791393
        w(6,2)=0.46791393
        w(6,3)=0.36076157
        w(6,4)=0.36076157
        w(6,5)=0.17132449
        w(6,6)=0.17132449
        w(7,1)=0.41795918
        w(7,2)=0.38183005
        w(7,3)=0.38183005
        w(7,4)=0.27970539
        w(7,5)=0.27970539
        w(7,6)=0.12948497
        w(7,7)=0.12948497
        w(8,1)=0.36268378
        w(8,2)=0.36268378
        w(8,3)=0.31370665
        w(8,4)=0.31370665
        w(8,5)=0.22238103
        w(8,6)=0.22238103
        w(8,7)=0.10122854
        w(8,8)=0.10122854
        w(9,1)=0.33023936
        w(9,2)=0.31234708
        w(9,3)=0.31234708
        w(9,4)=0.26061070
        w(9,5)=0.26061070
        w(9,6)=0.18064816
        w(9,7)=0.18064816
        w(9,8)=0.08127439
        w(9,9)=0.08127439
        w(10,1)=0.29552422
        w(10,2)=0.29552422
        w(10,3)=0.26926672
        w(10,4)=0.26926672
        w(10,5)=0.21908636
        w(10,6)=0.21908636
        w(10,7)=0.14945135
        w(10,8)=0.14945135
        w(10,9)=0.06667134
        w(10,10)=0.06667134

      End SUBROUTINE tw_vals