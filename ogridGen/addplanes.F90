!
!     Program that copies L-planes of single airfoil
!
      PROGRAM ADDPLANES
!
!
      INTEGER :: NJ,NK,NL,NLN
      INTEGER :: NGRID
      REAL(KIND=8) :: DS,TE,SP,YP
      REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: X,Y,Z,XO,YO,ZO 
      CHARACTER(80) :: OFILE,IFILE
      LOGICAL :: IOERR
!
!     Pulling this file straight from the hypgen output
      IOERR = .FALSE.
      OPEN(1,file='plot3d.dat',status='OLD',form='UNFORMATTED',err=10)
      READ(1,end=10,err=10) NGRID
      READ(1,end=10,err=10) NJ,NK,NL

      ALLOCATE( X(NJ,NK,NL), Y(NJ,NK,NL), Z(NJ,NK,NL) ) 
      READ(1,end=10,err=10) X,Y,Z
!
!     We are swapping out K and L 
      NLN = 3
      ALLOCATE( XO(NJ,NL,NLN), YO(NJ,NL,NLN), ZO(NJ,NL,NLN) )

      DO J=1,NJ
        DO K=1,NL
          DO L=1,NLN  
            XO(J,K,L) = X(J,1,K)
            YO(J,K,L) = Z(J,1,K)
            ZO(J,K,L) = 2-L
          ENDDO
        ENDDO 
      ENDDO 

      OPEN(2,file='grid.in',status='REPLACE',form='UNFORMATTED') 
      WRITE(2) 1 
      WRITE(2) NJ,NL,NLN
!     Switch z and y to 'rotate' 
      WRITE(2) XO,ZO,YO

      CLOSE(2) 
      CLOSE(1)
 10   CONTINUE
      END
