  SUBROUTINE INFLOW(IV,UX,UY,NNODE,VNPNT,NELEM_PP,NPOIN,INTMA,&
     &         CINF,NBOUN,NBNOI,BSIDO,DISNF,ITIME,PS,R,M,&
     &       VSPACE_FIRST,VSPACE_LAST)
! 
! *** THIS ROUTINE RESETS ALL DISCONTINUOUS NODES AT AN INFLOW BOUNDARY TO THE
! *** PRESCRIBED CONDITIONS AT INFINITY (FROM CINF)
!
      IMPLICIT NONE 
! 
      INTEGER IV,NNODE,VNPNT,NELEM_PP,NPOIN,NBOUN,NBNOI,IP1,IP2,TEST,IPT
      INTEGER IB,IE,IN,ITIME
      INTEGER INTMA(NNODE,NELEM_PP),BSIDO(NBNOI,NBOUN) 
      REAL nfw,Pw,Tw,RHOw,BETAw,C0w,C1w,C2w 
      REAL UX,UY,Tin,Pin,U0,V0,R,NA,PI,M,RHO,BETA,C0,C1,C2,UDASH,VDASH,nf 
      REAL CINF(4),PS(NPOIN)
      INTEGER VSPACE_FIRST,VSPACE_LAST
      REAL DISNF(3,VSPACE_FIRST:VSPACE_LAST,NELEM_PP) 

!
      PARAMETER (NA=6.022E+26,PI=3.1416)
!
! *** CONSTRUCT THE MAXWELLIAN nf FOR THIS POINT IN VELOCITY SPACE
! *** BASED ON THE PRESCRIBED INFLOW CONDITIONS
!
        Tin=CINF(1)
        Pin=CINF(2)*(10**5)
        U0=CINF(3)
        V0=CINF(4)
        RHO=Pin/(R*Tin)
        BETA=SQRT(RHO/(2*Pin))
        C0=(RHO*NA)/M
        C1=(BETA**2)/PI
        UDASH=UX-U0
        VDASH=UY-V0
        C2=(UDASH**2+VDASH**2)
        nf=C0*C1*EXP(-(BETA**2)*C2)
!
! *** LOOP OVER ALL BOUNDARY SIDES
!
       DO 1001 IB=1,NBOUN
        IP1=BSIDO(1,IB)
        IP2=BSIDO(2,IB)
        TEST=BSIDO(4,IB)
       IF(TEST.EQ.1)THEN     !INFLOW BOUNDARY
! *** FIND ALL THE DISCONTINUOUS NODES SURROUNDING THIS MESH NODE
       DO 1002 IE=1,NELEM_PP
       DO 1003 IN=1,NNODE
          IPT=INTMA(IN,IE)
       IF((IPT.EQ.IP1).OR.(IPT.EQ.IP2))THEN
            DISNF(IN,IV,IE)=nf
       ENDIF
!
 1003 CONTINUE
 1002 CONTINUE
       ENDIF
!
 1001 CONTINUE
! 
          RETURN 
! 
	  END 
