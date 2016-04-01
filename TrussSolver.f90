! Christopher Wenger
! Aere 361 Truss Solver Code
!This program will calculate the unknown displacement and forces of the truss.
!It accepts input from a text or dat file and intilizes the truss and begins to solve for 
!the unknown values

!                FILE INPUT FORMAT                               
!	 Total Number of pins  						    
! 	 Total Number of members							
!	 Pins the i-th member is contected to 
!			|
!	 Starting pin 	|    Ending Pin     
!			|
! 	 Until the n-th member is contected to 
! 	 i-th Pin's locations 
!			|
!		x 	|	 y          
!	 Until the n-th member
!	 i-th Member's Cross-section Area , Young's Modulus 
!			|
!			|
!	 Until the n-th member
!	 Known displacements and forces 
!			|
!	 0 --> given displacement  value of displacement
!	 1 --> given force 	   value of force
!	 Until 2 * total number of pins
!		  
!		CALL ORDER
!	TrussDef --> defines the truss 
!	MemElePropes --> defines the member and element properties
!	mechanics --> defines the stiffness matrix
!	flagger --> swaps the unknown and know elements in the matrices
!	DECOMP
!	SOLVE
!	PostProcessor --> organizes data and calculates stress, strain and load
!       dataProcessor --> Prints data to dat file
!       graphicsProcessor --> Calculates new x and y points and puts them in dat file for graphic software
PROGRAM TrussSolver
 IMPLICIT NONE
 INTEGER :: NMEM,NPIN,NDIM=0,N=0,ierror = 0,i,j,counter,low = 1
 REAL(KIND=8) :: COND
 REAL(KIND=8),ALLOCATABLE ::x(:),y(:),xNew(:),yNew(:),length(:),c(:),s(:),area(:),youngs(:)
 REAL(KIND=8),ALLOCATABLE :: given(:),v(:),B(:),A(:,:),u(:),f(:),k(:),stiff(:,:),WORK(:),strain(:),stress(:),Load(:)
 INTEGER, ALLOCATABLE :: IPVT(:),NCON(:,:)
 CHARACTER(20) :: filename
 WRITE(*,*)"Enter in file name with declaration"
 READ(*,*)filename
 OPEN(Unit=1,File=filename,Status='OLD',Action='READ',Iostat=ierror)
 READ(1,*) NPIN
 READ(1,*) NMEM
 ALLOCATE(NCON(low:NMEM,2),x(low:NPIN),y(low:NPIN),xNew(low:NPIN),yNew(low:NPIN))
 ALLOCATE(length(low:NMEM),c(low:NMEM),s(low:NMEM),area(low:NMEM),youngs(low:NMEM))
 ALLOCATE(given(low:2*NPIN),v(low:2*NPIN),u(low:2*NPIN),f(low:2*NPIN),k(low:NMEM),stiff(low:2*NPIN,low:2*NPIN))
 ALLOCATE(b(low:2*NPIN),A(low:2*NPIN,low:2*NPIN),IPVT(low:2*NPIN),WORK(low:2*NPIN),strain(low:NMEM),stress(low:NMEM),Load(low:NMEM))
 NDIM = 2*NPIN
 N = NDIM !Sets N to the size of NDIM
 !Defines the truss
 CALL TrussDef(NCON,x,y,NMEM,NPIN)
 !Gets the member and element properties
 CALL MemEleProps(NCON,NMEM,NPIN,x,y,length,c,s,area,youngs,v,given)
 !Defines the stiffness matrix of the truss
 CALL mechanics(NCON,c,s,NMEM,NPIN,length,area,youngs,k,stiff)
 !Lets organizes the matrices
 CALL flagger(NPIN,given,stiff,v,B,A)
 CALL DECOMP(NDIM,N,COND,IPVT,WORK,A)
	   IF (COND .NE. 1.0E32) THEN
	   CALL SOLVE (NDIM,N,B,IPVT,A)
	   ELSE
		WRITE(*,*)"ERROR THERE BOB"
	   END IF
 CALL PostProcessor(NPIN,NMEM,NCON,given,B,v,u,f,c,s,area,youngs,length,strain,stress,Load)
 CALL dataProcessor(NMEM,NPIN,u,f,stress,strain,Load)
 CALL graphicsProcesor(NCON,NMEM,NPIN,x,y,u,xNew,yNew)
 CLOSE(1)
 END PROGRAM
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Truss Definition!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE TrussDef(NCON,x,y,NMEM,NPIN)
  IMPLICIT NONE
 INTEGER,INTENT(IN) :: NMEM,NPIN
 REAL(KIND=8), INTENT(OUT) :: x(NPIN),y(NPIN)
 INTEGER, INTENT(OUT) :: NCON(NMEM,2)
 INTEGER :: i
!collects the pins the members are connected to
 DO i = 1,NMEM
 READ(1,*)NCON(i,1:2)
 END DO
!collects the pins' location (x,y)
 DO i = 1,NPIN
 READ(1,*)x(i),y(i)
 END DO
 END SUBROUTINE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Member/Element Properties!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE MemEleProps(NCON,NMEM,NPIN,x,y,length,c,s,area,youngs,v,given)
  IMPLICIT NONE
 INTEGER, INTENT(IN) :: NMEM, NPIN, NCON(NMEM,2)
 REAL(KIND=8), INTENT(IN) :: x(NPIN),y(NPIN)
 REAL(KIND=8), INTENT(OUT) :: length(NMEM),c(NMEM),s(NMEM), area(NMEM),youngs(NMEM)
 REAL(KIND=8), INTENT(OUT) :: v(2*NPIN),given(2*NPIN)
 REAL(KIND=8) Dx,Dy
 INTEGER :: i
 !collects the change in position for each member
 !then calculates the length and sin and cos of the member
 DO i = 1,NMEM
 Dx = x(NCON(i,2))-x(NCON(i,1))
 Dy = y(NCON(i,2))-y(NCON(i,1))
 length(i) = SQRT(Dx*Dx+Dy*Dy)
 c(i) = Dx/length(i)
 s(i) = Dy/length(i)
 END DO 
 !collects the area and young's modulus of each member
 DO i = 1,NMEM
 READ(1,*)area(i),youngs(i)
 END DO
 !collects the given forces and displacements at the pins
 DO i = 1, 2*NPIN
 READ(1,*)given(i),v(i)
 END DO
 END SUBROUTINE

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Flagging some variables !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE flagger(NPIN,given,stiff,v,B,A)
   IMPLICIT NONE
  INTEGER, INTENT(IN) ::NPIN
  REAL(KIND=8), INTENT(IN) :: given(2*NPIN), v(2*NPIN),stiff(2*NPIN,2*NPIN)
  REAL(KIND=8), INTENT(OUT) :: B(2*NPIN),A(2*NPIN,2*NPIN)
  INTEGER :: i,j
  REAL(KIND=8) :: scal = 100000.0
 DO i = 1,2*NPIN
	B(i) = 0.0 
 	DO j = 1,2*NPIN
	A(i,j) = stiff(i,j)
	END DO
 END DO
  DO j = 1,2*NPIN
	IF(given(j) .EQ. 0) THEN
	   DO i = 1,2*NPIN
		B(i) = B(i) - stiff(i,j)*v(j)
	   END DO
	   DO i = 1,2*NPIN
   		A(i,j) = 0.0
	   END DO
	   A(j,j) = -1.0
	ELSE
	   B(j) = B(j) +v(j)/scal
	END IF
  END DO
 END SUBROUTINE
 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Mechanics!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 SUBROUTINE mechanics(NCON,c,s,NMEM,NPIN,length,area,youngs,k,stiff)
  IMPLICIT NONE
 INTEGER, INTENT(IN) :: NMEM, NPIN, NCON(NMEM,2)
 INTEGER :: i,j,b,e,tbm1,tb,tem1,te
 REAL(KIND=8), INTENT(IN) :: c(NMEM),s(NMEM),length(NMEM), area(NMEM),youngs(NMEM)
 REAL(KIND=8), INTENT(OUT) :: k(NMEM),stiff(2*NPIN,2*NPIN)
 REAL(KIND=8) :: kcs, kc2, ks2
DO i = 1, 2*NPIN
DO j = 1, 2*NPIN
	stiff(i,j) = 0.0
END DO
END DO
!Lets calculate the stiffness matrix
	DO i =1,NMEM
		b = NCON(i,1)
		e = NCON(i,2)
		tbm1 = 2*b - 1
		tb = tbm1 + 1
		tem1 = 2*e - 1
		te = tem1 + 1
		k(i) = ((area(i)*youngs(i))/length(i))/10E4
	
		kcs = k(i)*c(i)*s(i)
		kc2 = k(i)*c(i)*c(i)
		ks2 = k(i)*s(i)*s(i)
		stiff(tbm1,tbm1) = stiff(tbm1,tbm1) + kc2
		stiff(tbm1,tb) = stiff(tbm1,tb) + kcs
		stiff(tbm1,tem1) = stiff(tbm1,tem1) - kc2
		stiff(tbm1,te) = stiff(tbm1,te) - kcs
		stiff(tb,tbm1) = stiff(tb,tbm1) + kcs
		stiff(tb,tb) = stiff(tb,tb) + ks2
		stiff(tb,tem1) = stiff(tb,tem1) - kcs
		stiff(tb,te) = stiff(tb,te) -ks2
		stiff(tem1,tbm1) = stiff(tem1,tbm1) -kc2
		stiff(tem1,tb) = stiff(tem1,tb) - kcs
		stiff(tem1,tem1) = stiff(tem1,tem1) + kc2
		stiff(tem1,te) = stiff(tem1,te) +kcs
		stiff(te,tbm1) = stiff(te,tbm1) - kcs
		stiff(te,tb) = stiff(te,tb) - ks2
		stiff(te,tem1) = stiff(te,tem1) + kcs
		stiff(te,te) = stiff(te,te) + ks2
	END DO
 END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Post Processor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PostProcessor(NPIN,NMEM,NCON,given,B,v,u,f,c,s,area,youngs,length,strain,stress,Load)
IMPLICIT NONE
 INTEGER, INTENT(IN) :: NPIN,NMEM,NCON(NMEM,2)
 REAL(KIND=8), INTENT(IN) :: given(2*NPIN), v(2*NPIN), B(2*NPIN),c(NMEM),s(NMEM),youngs(NMEM),length(NMEM),area(NMEM)
 REAL(KIND=8), INTENT(OUT) :: u(2*NPIN),f(2*NPIN),strain(NMEM),stress(NMEM),Load(NMEM)
 INTEGER :: i,b2,e,tbm1,tb,tem1,te
 REAL(KIND=8) :: scal = 100000
 	DO i =1,2*NPIN
	   IF (given(i) .EQ. 0) THEN
		u(i) = v(i)
	        f(i) = B(i)
           ELSE
		u(i) = B(i)
 		f(i) = v(i)/scal
	   END IF
	END DO
	DO i =1,NMEM
		b2 = NCON(i,1)
		e = NCON(i,2)
		tbm1 = 2*b2 - 1
		tb = tbm1 + 1
		tem1 = 2*e - 1
		te = tem1 + 1
		strain(i) = 1/length(i)*(c(i)*(u(tem1)-u(tbm1))+s(i)*(u(te)-u(tb)))
		stress(i) = youngs(i)*strain(i)
		Load(i) = area(i)*stress(i)
	END DO
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Data Processor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE dataProcessor(NMEM,NPIN,u,f,stress,strain,Load)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NMEM, NPIN
REAL(KIND=8), INTENT(IN) :: u(2*NPIN), f(2*NPIN), stress(NMEM), strain(NMEM), Load(NMEM)
REAL(KIND=8) :: scal = 100000.0
INTEGER :: i,counter,ierror = 0
OPEN(Unit=2,File="Output.dat",Status='OLD',Action='WRITE',Iostat=ierror)
 counter = 2
 DO i = 1,NMEM
    WRITE(2,53) "For Member", i
    WRITE(2,56) "Strain = ", strain(i)
    WRITE(2,57) "Stress = ", stress(i)
    WRITE(2,57) "Load =", Load(i)
 END DO
 DO i = 1,NPIN
    WRITE(2,53) "For pin", i
    WRITE(2,55)"Displacements = ", u(counter-1), u(counter)
    WRITE(2,54)"Forces  = ", f(counter-1)*scal, f(counter)*scal
	counter = counter + 2
 END DO
 53 FORMAT(A15,I2)
 54 FORMAT(A15,F10.3,X,F10.3)
 55 FORMAT(A15,F10.6,X,F10.6)
 56 FORMAT(A15,F10.6)
 57 FORMAT(A15,F10.3)
 CLOSE(2)
END SUBROUTINE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Graphic Processor!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE graphicsProcesor(NCON,NMEM,NPIN,x,y,u,xNew,yNew)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NPIN,NMEM,NCON(NMEM,2)
REAL(KIND=8), INTENT(IN) :: x(NPIN),y(NPIN), u(2*NPIN)
REAL(KIND=8), INTENT(OUT) :: xNew(NPIN),yNew(NPIN)
INTEGER :: i,counter,ierror = 0
OPEN(Unit=3,File="Plot.dat",Status='OLD',Action='WRITE',Iostat=ierror)
 DO i = 1,NPIN
    WRITE(3,58) x(i),y(i)
 END DO
 counter = 2
 DO i = 1,NPIN
    xNew(i) = x(i) + u(counter-1)
    yNew(i) = y(i) + u(counter)
	counter = counter + 2
    WRITE(3,58) xNew(i), yNew(i)
 END DO
 DO i = 1,NMEM
    WRITE(3,59) NCON(i,1),NCON(i,2)
 END DO
  58 FORMAT(F10.3,X,F10.3)
  59 FORMAT(I2,X,I2)
 CLOSE(3)
END SUBROUTINE
