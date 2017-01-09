!
!     Program that adds points to the finite trailing edge
      PROGRAM ADDENDPOINTS
!
!
      INTEGER :: NJ,NK,NL,NJN,NKN,NLN,NPTS
      INTEGER :: NGRID
      REAL(KIND=8) :: DS,TE,SP,YP
      REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: X,Y,Z,XO,YO,ZO 
      CHARACTER(80) :: OFILE,IFILE
      LOGICAL :: IOERR
!
!     We know there is only one grid, pulling this straight from
!     airfoil bezier code
!
      IOERR = .FALSE.
      OPEN(1,file='airfoil.p3d',status='OLD',form='UNFORMATTED',err=10)
      READ(1,end=10,err=10) NGRID
      READ(1,end=10,err=10) NJ,NK,NL

      ALLOCATE( X(NJ,NK,NL), Y(NJ,NK,NL), Z(NJ,NK,NL) ) 
      READ(1,end=10,err=10) X,Y,Z
!
!     Get trailing edge spacing 
!
      DS = SQRT( (X(2,1,1)-X(1,1,1))**2 + (Y(2,1,1)-Y(1,1,1))**2 )
      TE = SQRT( (X(NJ,1,1)-X(1,1,1))**2 + (Y(NJ,1,1)-Y(1,1,1))**2 )

      NPTS = TE/(DS*0.9D0)    
      SP = TE/NPTS
      NJN = NJ + NPTS

      ALLOCATE( XO(NJN,1,1), YO(NJN,1,1), ZO(NJN,1,1) )

      DO J = 1,NJ
        XO(J,1,1) = X(J,1,1)
        YO(J,1,1) = Y(J,1,1)
        ZO(J,1,1) = Z(J,1,1)
      ENDDO 

!     Add points at the end 
      DO I=1,NPTS
        YO(NJ+I,1,1) = YO(NJ+I-1,1,1) - SP
        XO(NJ+I,1,1) = X(NJ,1,1)
        ZO(NJ+I,1,1) = Z(NJ,1,1)
      ENDDO 

      OPEN(2,file='surf.i',status='REPLACE',form='UNFORMATTED') 
      WRITE(2) 1 
      WRITE(2) NJN,NK,NL
!     Switch z and y to 'rotate' 
      WRITE(2) XO,ZO,YO

      CLOSE(2) 
      CLOSE(1)
 10   CONTINUE
      END
