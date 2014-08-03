      SUBROUTINE IMPMM(I)
*
*
*       Multiple collision or merger search.
*       ------------------------------------
*
      INCLUDE 'common6.h'
!      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
!     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
!     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
!      COMMON/CLUMP/   BODYS(NCMAX,5),T0S(5),TS(5),STEPS(5),RMAXS(5),
!     &                NAMES(NCMAX,5),ISYS(5)
!      CHARACTER*8  WHICH1
      REAL*8  XX(3,3),VV(3,3)
 !     INTEGER LISTQ(100)
 !     SAVE LISTQ,QCHECK
 !     DATA IZARE,LISTQ(1),QCHECK /0,0,0.0D0/
*
*
*       Set index of KS pair & both components of c.m. body #I.
      IPAIR = I - N
      I1 = 2*IPAIR - 1
      I2 = I1 + 1
      !NTTRY = NTTRY + 1
      JCOMP = IFIRST
      !KS2 = 0
      RMAX2 = 1.0
      TTOT = TIME + TOFF
      RI2 = (X(1,I) - RDENS(1))**2 + (X(2,I) - RDENS(2))**2 +
     &                               (X(3,I) - RDENS(3))**2
*
*       Search c.m. neighbours if binary has at most two perturbers.
      J1 = I1
      IF (LIST(1,J1).LE.2) J1 = I
      NNB2 = LIST(1,J1) + 1
*
*       Find the dominant body (JCOMP) and nearest perturber (JMAX).
      RCRIT2 = 1.0D+04*RMIN2
      ITER = 0
    5 PERT1 = 0.0
      PERT2 = 0.0
      NP = 0
      DO 10 L = 2,NNB2
          J = LIST(L,J1)
          RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                  (X(3,I) - X(3,J))**2
          IF (J1.EQ.I.AND.RIJ2.GT.RCRIT2) GO TO 10
          NP = NP + 1
          JLIST(NP) = J
          PERT = BODY(J)/(RIJ2*SQRT(RIJ2))
          IF (PERT.GT.PERT2) THEN 
              IF (PERT.GT.PERT1) THEN
                  RJMIN2 = RIJ2
                  JMAX = JCOMP
                  JCOMP = J
                  PERT2 = PERT1
                  PERT1 = PERT
              ELSE
                  JMAX = J
                  PERT2 = PERT
                  RMAX2 = RIJ2
              END IF
          END IF
   10 CONTINUE
*
*       Perform an iteration to ensure at least two perturbers (< 10 times).
      IF (NP.LT.2) THEN
          RCRIT2 = 4.0*RCRIT2
          ITER = ITER + 1
          IF (ITER.LT.10) GO TO 5
      END IF
*
      RDOT = (X(1,I) - X(1,JCOMP))*(XDOT(1,I) - XDOT(1,JCOMP)) +
     &       (X(2,I) - X(2,JCOMP))*(XDOT(2,I) - XDOT(2,JCOMP)) +
     &       (X(3,I) - X(3,JCOMP))*(XDOT(3,I) - XDOT(3,JCOMP))
*
*
*       Only accept inward motion 
      IF (RDOT.GT. 0.0) GO TO 100
*
*       Include impact parameter test to distinguish different cases.
      A2 = (XDOT(1,I) - XDOT(1,JCOMP))**2 + 
     &     (XDOT(2,I) - XDOT(2,JCOMP))**2 +
     &     (XDOT(3,I) - XDOT(3,JCOMP))**2
      RIJ = SQRT(RJMIN2)
      A3 = 2.0/RIJ - A2/(BODY(I) + BODY(JCOMP))
      SEMI1 = 1.0/A3
      A4 = RDOT**2/(SEMI1*(BODY(I) + BODY(JCOMP)))
      ECC1 = SQRT((1.0D0 - RIJ/SEMI1)**2 + A4)
      PMIN = SEMI1*(1.0D0 - ECC1)
*
*       Set semi-major axis, eccentricity & apocentre of inner binary.
      SEMI = -0.5D0*BODY(I)/H(IPAIR)
      A0 = SEMI
      ECC2 = (1.0D0 - R(IPAIR)/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
      ECC = SQRT(ECC2)
      APO = ABS(SEMI)*(1.0 + ECC)
      PERIOD = TWOPI*A0*SQRT(A0/(BODY(I1)+BODY(I2)))
*
*       Quit on hyperbolic orbit with large impact parameter.
      IF (ECC1.GT.1.0.AND.PMIN.GT.50.0*SEMI) GO TO 100

*       
*       MM modifying these diagnostics
      ! binding energy of inner binary
      !EB = BODY(I1)*BODY(I2)*H(IPAIR)/BODY(I)
*       Print optional diagnostics for strong three-body interactions.
      !IF (KZ(45).GT.3.AND.BODY(I).GT.2.0D-04.AND.PMIN.LT.3.0*SEMI) THEN
      EB = H(IPAIR)
      IF (KZ(45).GT.3 .AND. EB.LT.-ECLOSE) THEN
!      pick which body is primary and companion based on mass 
         IF (BODY(I1).GT.BODY(I2)) THEN
            IPRIM = I1
            ICOMP = I2
         ELSE
            IPRIM = I2
            ICOMP = I1
         END IF
         ! APPLY THE PERICENTER CRITERIA AND THAT PRIM=MBH
        IF(PMIN.LT.3.D0*A0 .AND. BODY(IPRIM).EQ.BODY1) THEN
            ! if this is a repeat entry 
            !  (ie less than 1 orb period) then skip
            IF((TIME+TOFF)-TWRITE.LT.PERIOD) GO TO 100
            ! vinf of interloper
            HI = 0.5*A2 - (BODY(I) + BODY(JCOMP))/RIJ
            VINF = 0.0
            IF (HI.GT.0.0) VINF = SQRT(2.0*HI)*VSTAR
            ! Now write to the file
            NWRITE = NWRITE + 1
            WRITE (50,*) 'IMPACT',NWRITE,
     &           (TIME+TOFF),GAMMA(IPAIR),
     &           NAME(IPRIM), BODY(IPRIM),KSTAR(IPRIM),RADIUS(IPRIM),
     &           NAME(ICOMP), BODY(ICOMP),KSTAR(ICOMP),RADIUS(ICOMP),
     &           A0,ECC,PERIOD,
     &           NAME(JCOMP),BODY(JCOMP),KSTAR(JCOMP),RADIUS(JCOMP),
     &           PMIN,ECC1,VINF
            CALL FLUSH(50)
            TWRITE = TIME+TOFF
            ! now write another file with the positions and velocities
            ! Copy coordinates and velocities to local variables.
            WRITE (51,*) 'IMPACT',NWRITE,(TIME+TOFF),
     &           BODY(IPRIM),RADIUS(IPRIM),X(:,IPRIM),XDOT(:,IPRIM),
     &           BODY(ICOMP),RADIUS(ICOMP),X(:,ICOMP),XDOT(:,ICOMP),
     &           BODY(JCOMP),RADIUS(JCOMP),X(:,JCOMP),XDOT(:,JCOMP)
            CALL FLUSH(51)
         END IF
      END IF


  100 RETURN !IF (IPHASE.NE.8) JCMAX = 0
*
      !RETURN
*
      END
