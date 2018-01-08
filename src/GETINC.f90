      SUBROUTINE GETINC(NNODE,MMAT,RHS,DELUN, NELEM,RESIDUAL,& 
     &         maxNELEM_PP,VNPNT,DISNF,UMEAN,NFO,MU,IV,DTE,&
     &         VSPACE_FIRST,VSPACE_LAST) 
! 
      INTEGER NELEM,NNODE,IV,VNPNT 
      REAL RHS(1,3,NELEM),DELUN(1,3,NELEM) 
      REAL MMAT(NNODE,maxNELEM_PP) 
      REAL RESIDUAL,RESIDUAL1,COLL,NF,DTE
      REAL UMEAN(NELEM),NFO(NNODE,NELEM) 
      REAL MU(NNODE,NELEM)
      INTEGER VSPACE_FIRST,VSPACE_LAST
      REAL DISNF(3,VSPACE_FIRST:VSPACE_LAST,NELEM_PP) 
! 

! 
! *** CALCULATES THE INCREMENTS DELUN 
! 
      RESIDUAL=0.0 
! 
! *** LOOP OVER EACH ELEMEMT 
! 
      DO 1000 IE=1,NELEM 
! 
! *** LOOP OVER EACH NODE WITHIN THE ELEMENT 
!      
      DO 2000 IN=1,NNODE      
! *** COMPUTE THE HALFTIMESTEP BGK COLLISION TERM
!
!     NF=DISNF(IN,IV,IE)-UMEAN(IE)
      NF=DISNF(IN,IV,IE)
      COLL=MU(IN,IE)*(NFO(IN,IE)-NF)
!	  FIND THE MAX INCREMENT 
      RESIDUAL1=MMAT(IN,IE)*RHS(1,IN,IE)+COLL*DTE 
      IF(ABS(RESIDUAL1).GT.ABS(RESIDUAL))RESIDUAL=RESIDUAL1 
      DELUN(1,IN,IE)=DELUN(1,IN,IE)+MMAT(IN,IE)*RHS(1,IN,IE)+COLL*DTE 
 2000 CONTINUE 
 1000 CONTINUE
! 
      RETURN 
      END 
