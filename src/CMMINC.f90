        SUBROUTINE CMMINC(NNODE ,CMMAT ,RHS ,DELUN, NELEM,RESIDUAL,& 
     &           maxNELEM_PP,VNPNT,DISNF,UMEAN,NFO,MU,IV,DTE,&
     &           VSPACE_FIRST,VSPACE_LAST) 
! 
      INTEGER NELEM,NNODE,IV,VNPNT 
      REAL RHS(1,3,NELEM),DELUN(1,3,NELEM) 
      REAL CMMAT(3,3,maxNELEM_PP) 
      REAL RESIDUAL,RESIDUAL1,COLL,NF,UMEAN(NELEM),DTE
      REAL NFO(NNODE,NELEM),MU(NNODE,NELEM)
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
! *** MULTIPLY RHS BY THE ELEMENT CONSISTENT MASS MATRIX 
!      
      DO 2010 IN=1,NNODE
!
! *** COMPUTE THE HALFTIMESTEP BGK COLLISION TERM
!
      NF=DISNF(IN,IV,IE)
      COLL=MU(IN,IE)*(NFO(IN,IE)-NF)
      DELUN(1,IN,IE)=DELUN(1,IN,IE)+COLL*DTE   
      DO 2000 JN=1,NNODE   
      RESIDUAL1=CMMAT(IN,JN,IE)*RHS(1,IN,IE)+COLL*DTE 
      IF(RESIDUAL1.GT.RESIDUAL)RESIDUAL=RESIDUAL1 
      DELUN(1,IN,IE)=DELUN(1,IN,IE)+CMMAT(IN,JN,IE)*RHS(1,IN,IE) 
 2000 CONTINUE 
 2010 CONTINUE 
 1000 CONTINUE 
! 
      RETURN 
      END 
