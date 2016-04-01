 !nbasis = no. of basis functions
 !ndata = no. of paired observations
 PROGRAM Curvefit
 REAL(KIND=8) :: temp!temporary variable for finding NDIM
 REAL(KIND=8),ALLOCATABLE :: x(:),y(:),f(:),a(:)
 INTEGER :: i,ndata=0,nbasis,ierror = 0  
 CHARACTER(10) :: filename 
!allows the user to input a file from the keyboard
 WRITE(*,*)"Enter in the total number of basis funcations"
 READ(*,*)nbasis
 WRITE(*,*)"Enter in file name with declaraction which contains your (x,y) pairs"
 READ(*,*)filename
 OPEN(Unit=1,File=filename,Status='OLD',Action='READ',Iostat=ierror)
	IF(ierror .EQ. 0) THEN
!gets the size of the matrix from the file
	   DO
	     READ(1,*,Iostat=ierror)temp
	     IF(ierror.NE.0)EXIT !Exits loop if there is an error in file	     
           ndata = ndata + 1
   	   END DO
         END IF
 REWIND(Unit=1)!resets the position of reading the file
 ALLOCATE(x(ndata),y(ndata),a(nbasis))!Allocate memory
	   DO I=1,ndata
             READ(1,*) x(i),y(i)
           END DO
 CALL ZMatrix(x,y,ndata,nbasis,a)

	DO i = 1,nbasis
 	  WRITE(*,55) "a(",i,") =", a(i)
	  55 FORMAT(A2,I1,A3,F10.6)
	END DO
 END PROGRAM

 SUBROUTINE ZMatrix(x,y,ndata,nbasis,a)
 INTEGER, INTENT(IN) :: ndata,nbasis
 INTEGER :: i,j,k, IPVT(nbasis)
 REAL(KIND = 8) :: ZM(ndata,nbasis),ZT(nbasis,ndata),ZZT(nbasis,nbasis),Zy(nbasis),c(ndata),z(nbasis),a(nbasis)
 REAL(KIND = 8) :: WORK(nbasis),COND
 REAL(KIND = 8),INTENT(IN) :: x(ndata),y(ndata)
 N = nbasis
	DO i = 1,ndata
	  CALL Basis(x(i),z,nbasis)
	  c(i) = y(i)
          DO j = 1,nbasis
            ZM(i,j) = z(j)
          END DO
        END DO

	DO i= 1,nbasis
          DO j = 1,ndata
            ZT (i,j) = ZM (j,i)
          END DO
	END DO

	DO i = 1,nbasis
          DO j = 1,nbasis
            ZZT (i,j) = 0.0
             DO k = 1 , ndata
               ZZT (i,j) = ZT (i,k) * ZM (k,j) + ZZT (i,j)
             END DO
          END DO
	END DO

	DO i = 1,nbasis
     	  Zy (i) = 0.0
          DO k = 1 , ndata
            Zy (i) = ZT (i,k) * c (k) + Zy (i)
          END DO
	END DO
        CALL DECOMP(nbasis,N,COND,IPVT,WORK,ZZT)
	   IF (COND .NE. 1.0E32) THEN
	    CALL SOLVE (nbasis,N,Zy,IPVT,ZZT)
	END IF
	DO i = 1,nbasis
	  a(i) =  Zy(i)
	END DO
 END SUBROUTINE
!Program for Basis Functions
 SUBROUTINE Basis(x,z,nbasis)
 REAL(KIND=8),INTENT(OUT) :: z(nbasis)
 REAL(KIND=8),INTENT(IN) :: x
 z (1) = x
 z (2) = EXP(x)       ! basis function #1
 z (3) = COS(x)     ! basis function #2
 END SUBROUTINE
