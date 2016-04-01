 ! Author: Christopher Wenger
 ! Paddle Method calculates the pressure distribution on an air-foil
 ! and the lift of an air-foil when it is place in an inviscid parallel
 ! flow approaching the air-foil with and angle of attack

 !***********User Input types*************!
 ! 1) NACA number input			  !
 ! 2) air-foil data from .dat or .txt file !
 !****************************************!
 PROGRAM PaddleMethod
 IMPLICIT NONE
 !************CONSTANT VALUES*************!
  REAL(KIND=8) :: pi = 4.0d0*DATAN(1.0d0)
  REAL(KIND=8) :: one_twopi
  INTEGER :: NACA,NPOINTS = 0,USERTYPE,i,ierror = 0 !ierror set to 0 for success
  REAL(KIND=8) :: temp,density,Patm,COND
  REAL(KIND=8) :: chord,maxChamber,maxChamberPosition,maxThickness,AOA,Vf,Cl
  REAL(KIND=8),ALLOCATABLE :: x(:),y(:),length(:),c(:),s(:),x_bar(:),y_bar(:)
  REAL(KIND=8),ALLOCATABLE :: A(:,:),B(:),Vt(:),P(:),WORK(:)
  INTEGER,ALLOCATABLE :: NCON(:,:),IPVT(:)
  CHARACTER(10) :: filename 
  one_twopi = 0.5d0*pi
OPEN(Unit=2,File='Output.dat',Status='NEW',Action='WRITE',Iostat=ierror)
!
! User interface system through keyboard
!
  WRITE(*,*) "Enter in the input type"
  WRITE(*,*) "1 for NACA # input"
  WRITE(*,*) "2 for air-foil data from .dat or .txt file"
  READ(*,*) USERTYPE
!
! Control system for user interface, ref. user input types
! 
 IF(USERTYPE.EQ.1) THEN
	 WRITE(*,*) "Enter in 4-digit NACA #"
 	 READ(*,*) NACA
 	 WRITE(*,*) "Enter # of elements"
 	 READ(*,*) NPOINTS
 	 WRITE(*,*) "Enter in angle of attack in degs"
 	 READ(*,*) AOA
 	 maxChamber = NACA /1000/100.0d0
 	 maxChamberPosition = MOD(NACA/100,10)/10.0d0
 	 maxThickness = MOD(NACA,100)/100.0d0
	 AOA = AOA*(pi/180.0) ! Overwrite AOA from degrees to Radians
 !Allows the user to input a data file and set up the (x,y) points
 ELSEIF(USERTYPE.EQ.2) THEN
     WRITE(*,*)"Enter in file name with declaration"
     READ(*,*)filename
     WRITE(*,*) "Enter in angle of attack in degrees"
     READ(*,*) AOA
     OPEN(Unit=1,File=filename,Status='OLD',Action='READ',Iostat=ierror)
	IF(ierror .EQ. 0) THEN
!gets the size of the matrix from the file
	   DO
	     READ(1,*,Iostat=ierror)temp
	     IF(ierror.NE.0)EXIT !Exits loop if there is an error in file	     
           NPOINTS = NPOINTS + 1
   	   END DO
 	END IF
END IF
!
!DEFAULT SETTINGS
! 
  chord = 1.0d0;     	
  Vf = 1.0d0 ! free stream Velocity
  density = 1.225 !density at sea level kg/m^3
  Patm = 101325 !pressure at sea level Pa
 IF(NPOINTS.LE.0) THEN
    NPOINTS = ABS(NPOINTS)
 END IF
 !Makes sure number of points is even for symmetric upper and lower surfaces 
 IF(MOD(NPOINTS,2).NE.0) THEN
    NPOINTS = NPOINTS + 1
 END IF

!*************************************************!
!Section begins the processing part of the code   !
!*************************************************!

!*******************************ALLOCATE MEMORY*******************************!
 ALLOCATE(x(1:NPOINTS),y(1:NPOINTS),length(1:NPOINTS),c(1:NPOINTS),s(1:NPOINTS),NCON(1:NPOINTS,1:2))
 ALLOCATE(x_bar(1:NPOINTS),y_bar(1:NPOINTS),A(1:NPOINTS+1,1:NPOINTS+1),B(1:NPOINTS+1),Vt(NPOINTS),P(NPOINTS))
 ALLOCATE(IPVT(1:NPOINTS+1),WORK(1:NPOINTS+1))
!*****************************************************************************!
 !Generates (x,y) points based on the NACA information 
 !otherwise read (x,y) points from data file or text file
 IF(USERTYPE.EQ.1) THEN
 CALL AirfoilGen(chord,maxChamber,maxChamberPosition,maxThickness,x,y,NPOINTS)
 ELSEIF(USERTYPE.EQ.2) THEN
	REWIND(Unit=1)!resets the position of reading the fill
	DO i=1,NPOINTS
            READ(1,*) x(i),y(i)
        END DO
 	
 END IF
 CALL AirfoilProps(NCON,NPOINTS,x,y,length,c,s,x_bar,y_bar)
 !AOA = alpha 
 CALL NoPenetration(NPOINTS,NCON,x, y, x_bar, y_bar,c,s,AOA,Vf,pi,A,B)
 CALL KuttaCondition(NPOINTS,NCON,x, y, x_bar, y_bar,c,s,AOA,Vf,pi,A,B)
 CALL Decomp(NPOINTS+1,NPOINTS+1,COND,IPVT,WORK,A)
 CALL Solve(NPOINTS+1,NPOINTS+1,B,IPVT,A)
 CALL PostProc(NPOINTS,NCON,x, y, x_bar, y_bar,c,s,length,AOA,Vf,pi,B,Vt,density,Patm,P,Cl)
 WRITE(*,*) "Alpha at ",AOA*(pi/180), " degrees: CL = ", Cl
 DO i = 1,NPOINTS
 WRITE(2,*) P(i),Vt(i)
 CLOSE(1)
 CLOSE(2)
 
 END PROGRAM 


!
!Subroutine Generates Air-foil (x,y) points based on NACA digits
! 
 SUBROUTINE AirfoilGen(chord,maxChamber,maxChamberPosition,maxThickness,x,y,NPOINTS)
 IMPLICIT NONE
 INTEGER,INTENT(IN) :: NPOINTS
 REAL(KIND=8),INTENT(IN) :: chord,maxChamber,maxChamberPosition,maxThickness
 REAL(KIND=8),INTENT(OUT) :: x(NPOINTS),y(NPOINTS)
 INTEGER :: k,NbyTwo
 REAL(KIND=8) :: xloc,del_x,xMin,xMax,yc,yt,dycdx,theta

xMin = 0.0d0
xMax = chord
NbyTwo = NPOINTS/2
del_x = (xMax-xMin)/REAL(NbyTwo)
!
!Generates the lower surface of the air-foil
!
	DO k = 1, NbyTwo
	 xloc = chord - (k-1)*del_x
	     IF(xMin .LE. xloc .AND. xloc .LT. maxChamberPosition*chord) THEN
	     yc = maxChamber/(maxChamberPosition*maxChamberPosition)*(2.0d0*maxChamberPosition*xloc &
		  -xloc*xloc)
	     dycdx = maxChamber/(maxChamberPosition*maxChamberPosition)*(2.0d0*maxChamberPosition &
		  -2.0d0*xloc)
	     ELSEIF(xloc .GE.maxChamberPosition*chord .AND. xloc .LE. chord) THEN
	     yc = maxChamber/((1.0d0-maxChamberPosition)*(1.0d0-maxChamberPosition)) &
		  *((1.0d0-2.0d0*maxChamberPosition)+2.0d0*maxChamberPosition*xloc-xloc*xloc)
	     dycdx = maxChamberPosition/((1.0d0-maxChamberPosition)*(1.0d0-maxChamberPosition)) &
        	     *((1.0d0-2.0d0*maxChamberPosition)+2.0d0*maxChamberPosition-xloc*2.0d0)
	     END IF
	     !Thickness Distribution
	     yt = maxThickness/0.20d0*(0.296375d0*SQRT(xloc)-0.12635d0*xloc-0.35195d0 &
		  *(xloc*xloc)+0.283775d0*(xloc*xloc*xloc)-0.10185d0*(xloc*xloc*xloc*xloc))
	     !Calculates the final x and y positions
	     theta = ATAN(dycdx)
	     x(k) = xloc+yt*SIN(theta)
	     y(k) = yc-yt*COS(theta)
	END DO
!
!Generates the upper surface of the air-foil
!
	DO k = NbyTwo, NPOINTS
	xloc = xMin + (k - NbyTwo)*del_x
	     IF(xMin .LE. xloc .AND. xloc .LT. maxChamberPosition*chord) THEN
	     yc = maxChamber/(maxChamberPosition*maxChamberPosition)*(2*maxChamberPosition*xloc &
		  -xloc*xloc)
	     dycdx = maxChamber/(maxChamberPosition*maxChamberPosition)*(2.0d0*maxChamberPosition &
		  -2.0d0*xloc)
	     ELSEIF(xloc .GE.maxChamberPosition*chord .AND. xloc .LE. chord) THEN
	     yc = maxChamber/((1.0d0-maxChamberPosition)*(1.0d0-maxChamberPosition)) &
		  *((1.0d0-2.0d0*maxChamberPosition)+2.0d0*maxChamberPosition*xloc-xloc*xloc)
	     dycdx = maxChamberPosition/((1.0d0-maxChamberPosition)*(1.0d0-maxChamberPosition)) &
        	     *((1.0d0-2.0d0*maxChamberPosition)+2.0d0*maxChamberPosition-xloc*2.0d0)
	     END IF
	     !Thickness Distribution
	     yt = maxThickness/0.20d0*(0.296375d0*SQRT(xloc)-0.12635d0*xloc-0.35195d0 &
		  *(xloc*xloc)+0.283775d0*(xloc*xloc*xloc)-0.10185d0*(xloc*xloc*xloc*xloc))
	     !Calculates the final x and y positions
	     theta = ATAN(dycdx)
	     x(k) = xloc-yt*SIN(theta)
	     y(k) = yc+yt*COS(theta)
	END DO
 END SUBROUTINE

 !
 !Defines the air-foil properties
 !
 SUBROUTINE AirfoilProps(NCON,NPOINTS,x,y,length,c,s,x_bar,y_bar)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NPOINTS
 INTEGER, INTENT(OUT) :: NCON(NPOINTS,2)
 REAL(KIND=8), INTENT(IN) :: x(NPOINTS),y(NPOINTS)
 REAL(KIND=8), INTENT(OUT) :: length(NPOINTS),c(NPOINTS),s(NPOINTS),x_bar(NPOINTS),y_bar(NPOINTS)
 REAL(KIND=8) Dx,Dy
 INTEGER :: i,j,b,e
 !Connects all of the points to make the segments
 DO i = 1,NPOINTS-1
    NCON(i,1) = i
    NCON(i,2) = i+1
 END DO
 !Makes sure the segments are enclosed
 NCON(NPOINTS,1) = NPOINTS
 NCON(NPOINTS,2) = 1
 !collects the change in position for each segments
 !then calculates the length and sin and cos of the segments
 !calculates the midpoint positions of the segments
 DO i = 1,NPOINTS
 b = NCON(i,1)
 e = NCON(i,2)
 Dx = x(e)-x(b)
 Dy = y(e)-y(b)
 length(i) = SQRT(Dx*Dx+Dy*Dy)
 c(i) = Dx/length(i)
 s(i) = Dy/length(i)
 x_bar(i) = 0.5d0*(x(b)+x(e))
 y_bar(i) = 0.5d0*(y(b)+y(e))
 END DO
 END SUBROUTINE
 !
 !No-Penetration Condition
 !
 SUBROUTINE NoPenetration(NPOINTS,NCON,x, y, x_bar, y_bar,c,s,AOA,Vf,pi,A,B)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NPOINTS,NCON(NPOINTS,2)
 REAL(KIND=8), INTENT(IN) :: x(NPOINTS),y(NPOINTS),x_bar(NPOINTS),y_bar(NPOINTS),c(NPOINTS),s(NPOINTS)
 REAL(KIND=8), INTENT(IN) :: pi,Vf,AOA
 REAL(KIND=8), INTENT(OUT) :: A(NPOINTS+1,NPOINTS+1),B(NPOINTS+1)
 INTEGER :: i,j
 REAL(KIND=8) :: dx1,dx2,dy1,dy2,r1,r2,rij,rijp1,OneOver2Pi,vp,up,CTIMTJ,STIMTJ,beta
 OneOver2Pi = 1.0d0/(2.0d0*pi)
 
 DO i = 1,NPOINTS
	A(i,NPOINTS+1) = 0.0d0 !Setting elements equal to zero
	DO j = 1,NPOINTS
		dx1 = x_bar(i) - x(NCON(j,1))
		dx2 = x_bar(i) - x(NCON(j,2))
		dy1 = y_bar(i) - y(NCON(j,1))
		dy2 = y_bar(i) - y(NCON(j,2))
		rij = DSQRT(dx1*dx1+dy1*dy1)
		rijP1 = DSQRT(dx2*dx2+dy2*dy2)
		IF(i.EQ.j) THEN
			beta = pi
		ELSE
			beta = DATAN2((dy2*dx1-dx2*dy1),(dy2*dy1+dx2*dx1))
		END IF
		up = OneOver2Pi * (DLOG(rijP1/rij))
        vp = OneOver2Pi * beta
        CTIMTJ = c(i)*c(j) + s(i)*s(j)
        STIMTJ = s(i)*c(j) - s(j)*c(i)
        A(i,j) = up*STIMTJ + vp*CTIMTJ
        A(i,NPOINTS+1) = A(i,NPOINTS+1) + up*CTIMTJ - vp*STIMTJ
	END DO
	B(i) = Vf*DCOS(AOA)*s(i) - Vf*DCOS(AOA)*c(i) 
 END DO
 END SUBROUTINE
!
!Kutta Condition
!
 SUBROUTINE KuttaCondition(NPOINTS,NCON,x, y, x_bar, y_bar,c,s,AOA,Vf,pi,A,B)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NPOINTS,NCON(NPOINTS,2)
 REAL(KIND=8), INTENT(IN) :: x(NPOINTS),y(NPOINTS),x_bar(NPOINTS),y_bar(NPOINTS),c(NPOINTS),s(NPOINTS)
 REAL(KIND=8), INTENT(IN) :: pi,Vf,AOA
 REAL(KIND=8), INTENT(OUT) :: A(NPOINTS+1,NPOINTS+1),B(NPOINTS+1)
 INTEGER :: j,k
 REAL(KIND=8) :: dx1,dx2,dy1,dy2,r1,r2,rij,rijp1,OneOver2Pi,vp,up,CTKMTJ,STKMTJ,beta
 OneOver2Pi = 1.0d0/(2.0d0*pi)
 
 A(NPOINTS+1,NPOINTS+1) = 0.0d0 !setting elements equal to zero
 B(NPOINTS+1) = 0.0d0
	DO j = 1,NPOINTS
	DO k = 1,NPOINTS,NPOINTS-1
		dx1 = x_bar(k) - x(NCON(j,1))
		dx2 = x_bar(k) - x(NCON(j,2))
		dy1 = y_bar(k) - y(NCON(j,1))
		dy2 = y_bar(k) - y(NCON(j,2))
		rij = DSQRT(dx1*dx1+dy1*dy1)
		rijP1 = DSQRT(dx2*dx2+dy2*dy2)
		IF(k.EQ.j) THEN
			beta = pi
		ELSE
			beta = DATAN2((dy2*dx1-dx2*dy1),(dy2*dy1+dx2*dx1))
		END IF
		up = OneOver2Pi * (DLOG(rijP1/rij))
        vp = OneOver2Pi * beta
        CTKMTJ = c(k)*c(j) + s(k)*s(j)
        STKMTJ = s(k)*c(j) - s(j)*c(k)
		A(NPOINTS+1,NPOINTS+1) = A(NPOINTS+1,NPOINTS+1) + (1.0d0/(2.0d0*pi))*beta*CTKMTJ +(1.0d0/(2.0d0*pi))*STKMTJ*DLOG(rijP1/rij)
	END DO
 END DO
 B(NPOINTS+1) = -Vf*((cos(AOA)*c(1)+sin(AOA)*s(1))+(cos(AOA)*c(NPOINTS)+sin(AOA)*s(NPOINTS))) 
 END SUBROUTINE
!
! Post processor
!
SUBROUTINE PostProc(NPOINTS,NCON,x, y, x_bar, y_bar,c,s,length,AOA,Vf,pi,B,Vt,density,Patm,P,Cl)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: NPOINTS,NCON(NPOINTS,2)
 REAL(KIND=8), INTENT(IN) :: x(NPOINTS),y(NPOINTS),x_bar(NPOINTS),y_bar(NPOINTS),c(NPOINTS),s(NPOINTS),length(NPOINTS)
 REAL(KIND=8), INTENT(IN) :: pi,Vf,AOA,density,Patm
 REAL(KIND=8), INTENT(IN) :: B(NPOINTS+1)
 REAL(KIND=8), INTENT(OUT) :: Vt(NPOINTS),P(NPOINTS),Cl
 INTEGER :: j,i
 REAL(KIND=8) :: dx1,dx2,dy1,dy2,r1,r2,rij,rijp1,OneOver2Pi,vp,up,CTIMTJ,STIMTJ,CTIMTA,beta
 OneOver2Pi = 1.0d0/(2.0d0*pi)
 Cl = 0.0d0
DO i = 1,NPOINTS
	DO j = 1,NPOINTS
		dx1 = x_bar(i) - x(NCON(j,1))
		dx2 = x_bar(i) - x(NCON(j,2))
		dy1 = y_bar(i) - y(NCON(j,1))
		dy2 = y_bar(i) - y(NCON(j,2))
		rij = DSQRT(dx1*dx1+dy1*dy1)
		rijP1 = DSQRT(dx2*dx2+dy2*dy2)
		IF(i.EQ.j) THEN
			beta = pi
		ELSE
			beta = DATAN2((dy2*dx1-dx2*dy1),(dy2*dy1+dx2*dx1))
		END IF
		up = (DLOG(rijP1/rij))
        vp = OneOver2Pi * beta
        CTIMTJ = c(i)*c(j) + s(i)*s(j)
        STIMTJ = s(i)*c(j) - s(j)*c(i)
	CTIMTA = c(i)*DCOS(AOA) + s(i)*DSIN(AOA)
	Vt(i) = Vt(i) + Vf*CTIMTA + B(i)*OneOver2Pi *(beta*STIMTJ-CTIMTJ*up)+B(NPOINTS+1)*OneOver2Pi*(STIMTJ*up+beta*CTIMTJ)
	P(i) = Patm + 0.5d0*Vf*Vf - 0.5d0*Vt(i)*Vt(i)
	END DO
	Cl = Cl + 2*B(NPOINTS+1)/Vf*length(i)
 END DO
END SUBROUTINE
