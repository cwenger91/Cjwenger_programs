!Christopher Wenger
!Aere 361
!Program Linear solver
!This program reads an A matrix NxN and a B matrix 1xN and
!solves for that unknown varibles in the system of equations
PROGRAM LINEAR_SOLVER
     
     REAL(KIND=8) :: COND,temp !temporary variable for finding NDIM
     REAL(KIND=8), ALLOCATABLE :: A(:,:), WORK(:),B(:)
     INTEGER, ALLOCATABLE :: IPVT(:)
     INTEGER :: NDIM=0,N=0,low=1,I,ierror = 0 !ierror set to 0 for success

!allows the user to input a file from the keyboard
     CHARACTER(10) :: filename 
     WRITE(*,*)"Enter in file name with declaraction"
     READ(*,*)filename
     OPEN(Unit=1,File=filename,Status='OLD',Action='READ',Iostat=ierror)
	IF(ierror .EQ. 0) THEN
!gets the size of the matrix from the file
	   DO
	     READ(1,*,Iostat=ierror)temp
	     IF(ierror.NE.0)EXIT !Exits loop if there is an error in file	     
           NDIM = NDIM + 1
   	   END DO
	   N = NDIM !Setting N = NDIM
	   ALLOCATE(A(low:NDIM,low:NDIM),WORK(low:NDIM),B(low:NDIM),IPVT(low:NDIM),STAT=ierror)!Allocate memory
	   IF(ierror.EQ.0)THEN
	     REWIND(Unit=1)!resets the position of reading the file
	   END IF
	   DO I=1,NDIM
             READ(1,*) A(i,1:NDIM),B(i)
           END DO

	   CALL DECOMP(NDIM,N,COND,IPVT,WORK,A)
	   IF (COND .NE. 1.0E32) THEN

	    CALL SOLVE (NDIM,N,B,IPVT,A)

! WRITING CONDITION NUMBER


	    WRITE(*,54) "condition number = ",COND
	54 FORMAT(A19,F10.6)
! WRITING SOLUTION STORED IN B
	
	    DO I=1,N
	    WRITE(*,55) "Value ", I, " = ",B(I)
	55 FORMAT(A6,I1,A3,F10.6)
	    END DO
	    ELSE
	    WRITE(*,*)"Singularity"
	    END IF
      END IF  !End of file check if statement
      CLOSE(1)!closes file opened by program
END PROGRAM LINEAR_SOLVER

