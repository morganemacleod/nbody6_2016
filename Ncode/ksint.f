      SUBROUTINE KSINT(I1)
*
*
*
*       Regularized integration.
*       ------------------------
*
*
      INCLUDE 'common6.h'
      COMMON/BINARY/  CM(4,MMAX),XREL(3,MMAX),VREL(3,MMAX),
     &                HM(MMAX),UM(4,MMAX),UMDOT(4,MMAX),TMDIS(MMAX),
     &                NAMEM(MMAX),NAMEG(MMAX),KSTARM(MMAX),IFLAGM(MMAX)
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               ZP(NTMAX),ES(NTMAX),CZ(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      COMMON/SLOW0/  RANGE,ISLOW(10)
      COMMON/GAMDOT/  DGAM
      REAL*8  UI(4),UIDOT(4),XI(6),VI(6),FP(6),FD(6)
      LOGICAL IQ
      SAVE ITERM
      REAL*8  DE,DH,TINSP
*
*
*       Set second component, pair index & c.m. index.
      I2 = I1 + 1
      IPAIR = KVEC(I1)
      I = N + IPAIR
*
*       Define perturber membership & inverse c.m. mass.
      NNB0 = LIST(1,I1)
      BODYIN = 1.0/BODY(I)
*
*       Check for further unperturbed motion.
      IF (NNB0.EQ.0.AND.H(IPAIR).LT.0.0) THEN
*       Include possible eccentricity modulation of hierarchical binary.
          IF (NAME(I).LT.0) THEN
              IM = 0
              DO 1 K = 1,NMERGE
                  IF (NAMEM(K).EQ.NAME(I)) IM = K
    1         CONTINUE
              IF (IM.GT.0.AND.TIME.GT.TMDIS(IM)) THEN
                  SEMI = -0.5*BODY(I)/H(IPAIR)
                  ECC2 = (1.0 - R(IPAIR)/SEMI)**2 +
     &                                   TDOT2(IPAIR)**2/(BODY(I)*SEMI)
                  RP = SEMI*(1.0 - SQRT(ECC2))*(1.0 - 2.0*GAMMA(IPAIR))
                  IF (RP.LT.R0(IPAIR)) THEN
                      WRITE (77,4)  NAME(I1), SQRT(ECC2), RP, R0(IPAIR)
    4                 FORMAT (' TMDIS TERM    NM E RP R0',I7,1P,3E10.2)
                      CALL FLUSH(77)
                      INSTAB = INSTAB + 1
                      GO TO 90
                  ELSE IF (KZ(27).EQ.2) THEN
                      CALL ECCMOD(I,ITERM)
*       Update time on termination to prevent looping.
                      IF (ITERM.GT.0) THEN
                          T0(I1) = TIME
                          GO TO 90
                      END IF
                  END IF
*       Check any inner and outer circularizing binary using NAMEM & NAMEG.
                  DO 3K = 1,NCHAOS
                      IF (NAMEC(K).EQ.NZERO - NAMEM(IM)) THEN
*       Update unperturbed binary if T - TOSC > 10 Myr (cf. IMPACT & DECIDE).
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) THEN
                              T0(I1) = TIME
                              GO TO 90
                          END IF
                      END IF
                      IF (NAMEC(K).EQ.NAMEG(IM)) THEN
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) THEN
                              T0(I1) = TIME
                              GO TO 90
                          END IF
                      END IF
    3             CONTINUE
              END IF
          END IF
          DT1 = STEP(I1)
          CALL UNPERT(IPAIR)
*
*       Check updating of unperturbed relativistic KS binary.
          IF (KZ(11).NE.0.AND.LIST(1,I1).EQ.0) THEN
              CALL BRAKE4(I1,I2,DT1)
              IF (IPHASE.LT.0) GO TO 100
          END IF
*
*       Try re-initialize chain WD/BH system after dormant KS (#11 only).
          IF (KZ(11).NE.0.AND.NCH.EQ.0.AND.LIST(1,I1).GT.0) THEN
              IF (MIN(KSTAR(I1),KSTAR(I2)).GE.10) THEN
                  SEMI = -0.5*BODY(I)/H(IPAIR)
                  IF (SEMI.GT.10.0*RMIN) GO TO 100
                  WRITE (6,222)  TIME+TOFF, NAME(JCLOSE), KSTAR(I1),
     &                           KSTAR(I2), LIST(1,I1), GAMMA(IPAIR),
     &                           SEMI, R(IPAIR)
  222             FORMAT (' ACTIVATE CHAIN    T NMJ K* NP G A R ',
     &                                        F9.1,I7,3I4,1P,3E10.2)
                  KSPAIR = IPAIR
*       Restore unperturbed motion from BRAKE4 (NP = 0 fixes some problem).
                  IF (GAMMA(IPAIR).LT.1.0D-10) THEN
                      JCLOSE = 0
                      LIST(1,I1) = 0
                  END IF
                  KS2 = 0
*       Include case of binary as dominant perturber.
                  IF (JCLOSE.GT.N) THEN
                      KS2 = JCLOSE - N
                      IF (KS2.GT.KSPAIR) KS2 = KS2 - 1
                      JCOMP = JCLOSE
                      JP = JCLOSE - N
                      WRITE (6,223)  KSPAIR, KS2, JCLOSE, GAMMA(JP)
  223                 FORMAT (' BINARY PERT    KSP KS2 JCLOSE GJP ',
     &                                         2I4,I7,1P,E10.2)
                  ELSE
*       Avoid JCOMP > N & JCLOSE < N for spurious CALL KSTERM in DELAY.
                      JCOMP = 0
                  END IF
                  IPHASE = 8
*       Save KS parameters until end of block-step (JCMAX=0: no extra pert).
                  CALL DELAY(1,KS2)
                  JCMAX = 0
              END IF
          END IF
          GO TO 100
      END IF
*
*       Perform KS prediction of U, UDOT & H.
      CALL KSPRED(IPAIR,I1,I,BODYIN,UI,UIDOT,XI,VI)
*
*       Obtain the perturbing force & derivative.
      CALL KSPERT(I1,NNB0,XI,VI,FP,FD)
*
*       Save old radial velocity & relative perturbation and set new GAMMA.
      RDOT = TDOT2(IPAIR)
*     G0 = GAMMA(IPAIR)
      GI = SQRT(FP(1)**2 + FP(2)**2 + FP(3)**2)*R(IPAIR)**2*BODYIN
      GAMMA(IPAIR) = GI
*
*       Apply the Hermite corrector.
      CALL KSCORR(IPAIR,UI,UIDOT,FP,FD,TD2,TDOT4,TDOT5,TDOT6)
*
*       Increase regularization time-step counter and update the time.
      NSTEPU = NSTEPU + 1
      T0(I1) = TIME
*       Check for early return during termination (called from KSTERM).
      IF (IPHASE.NE.0) GO TO 100
*
*       Define useful scalars.
      RI = R(IPAIR)
      HI = H(IPAIR)
*
*       Initialize termination indicator and check for large perturbation.
      IQ = .FALSE.
      IF (GI.GT.0.25) GO TO 2
      IF (GI.LT.0.03) THEN
          JCOMP = 0
          GO TO 20
      END IF
      CALL FLYBY(I,ITERM)
      IF (ITERM.EQ.0.AND.KSTAR(I).LT.0) THEN
          IQ = .TRUE.
      END IF
      IF (ITERM.EQ.1) THEN
*       Delay chain regularization search until near end of block-step.
          IF (TIME + STEP(I1).GT.TBLOCK) THEN
              CALL IMPACT(I)
              ! MMadd
              IF(KZ(45) .GE. 3 ) CALL IMPMM(I)
              IF (IPHASE.GT.0) GO TO 100
          END IF
      ELSE IF (ITERM.EQ.2) THEN
          IQ = .TRUE.
          GO TO 20
      ELSE
          GO TO 20
      END IF
*
*       Find the dominant body for large perturbations.
    2 S = 4.0*STEP(I)
      FMAX = BODY(I)/RI**2
*       Initialize JCOMP for prediction and optional diagnostics in KSTERM.
      JCOMP = 0
      DO 10 L = 2,NNB0+1
          J = LIST(L,I1)
*       Only search bodies within twice the c.m. time-step.
          IF (STEP(J).GT.S) GO TO 10
*       Compare strong perturber and either component with current pair.
          DO 5 K = I1,I2
              RIJ2 = (X(1,J) - X(1,K))**2 + (X(2,J) - X(2,K))**2 +
     &                                      (X(3,J) - X(3,K))**2
              IF (BODY(J) + BODY(K).GT.RIJ2*FMAX) JCOMP = J
    5     CONTINUE
   10 CONTINUE
*
*       Set termination if strong perturber <= N forms dominant pair.
      IF (JCOMP.GT.0.OR.GI.GT.1.0) THEN
*       Check optional binary search.
*         IF (KZ(4).GT.0) THEN
*             DGAM = GI - G0
*             K = KZ(4)
*             CALL EVOLVE(IPAIR,K)
*         END IF
          IF (JCOMP.LE.N.OR.GI.GT.2.0) IQ = .TRUE.
      END IF
*
*       Check termination of hyperbolic encounter (R > R0 or R > 2*RMIN).
   20 IF (HI.GT.0.0D0.AND.NAME(I).GT.0) THEN
          IF ((RI.GT.R0(IPAIR).AND.GI.GT.GMAX).OR.RI.GT.2.0*RMIN.OR.
     &        (GI.GT.0.5.AND.TD2.GT.0.0)) THEN
*       Skip termination delay in case of velocity kick (cf. routine KSTERM).
              IF (HI.LT.100.0.OR.GI.GT.0.1.OR.RI.GT.5.0*RMIN) THEN
                  IF (TD2.GT.0.0) IQ = .TRUE.
              END IF
          END IF
      END IF
*
*       Choose basic regularized step using binding energy or lower limit.
      IF (ABS(HI).GT.ECLOSE) THEN
          W1 = 0.5/ABS(HI)
      ELSE
          W1 = R0(IPAIR)*BODYIN
          W2 = 0.5/ABS(HI)
          W1 = MIN(W1,W2)
          IF (RI.GT.R0(IPAIR)) W1 = W1*R0(IPAIR)/RI
*       Maximum square step for soft binaries & weak hyperbolic pairs.
          IF (HI.LT.0.0D0) THEN
*       Include case of hard binary with massive components or merger.
              W2 = -0.5/HI
              IF (NAME(I).LT.0) THEN
                  W1 = RI*BODYIN
              END IF
              W1 = MIN(W2,W1)
          END IF
      END IF
*
*       Include perturbation factor in predicted step.
      IF (GI.LT.0.0005) THEN
*       Use second-order expansion of cube root for small perturbations.
          W3 = 333.3*GI
          W2 = SQRT(W1)/(1.0 + W3*(1.0 - W3))
      ELSE
          W3 = 1.0 + 1000.0*GI
          W2 = SQRT(W1)/W3**0.3333
      END IF
*
*       Form new regularized step (include inertial factor).
      DTU = ETAU*W2
      DTU = MIN(1.2*DTAU(IPAIR),DTU)
*
*       Include convergence criterion DH = H'*DTU + H''*DTU**2/2 = 0.001*|H|.
      IF (GI.GT.1.0D-04) THEN
          DH = 1.0E-03*MAX(ABS(H(IPAIR)),ECLOSE)
          XF = 2.0*DH/ABS(HDOT2(IPAIR))
          YF = HDOT(IPAIR)/HDOT2(IPAIR)
          DTU1 = SQRT(XF + YF**2) - ABS(YF)
          DTU = MIN(DTU1,DTU)
      END IF
*
*       Check pericentre step reduction for perturbed spiral.
      IF (KSTAR(I).EQ.-2.AND.TD2.LT.0.0) THEN
          SEMI = -0.5*BODY(I)/HI
          IF (RI.LT.SEMI) THEN
*       Predict radial velocity for step DTU (note scaled coefficients).
              RD = ((ONE3*TDOT5*DTU + TDOT4)*DTU +
     &                                TDOT3(IPAIR))*DTU + 2.0*TD2
*       Adopt pericentre step with 1% safety factor (small TDOT4 is OK).
              IF (RD.GT.0.0) THEN
                  A2 = 0.5*TDOT3(IPAIR)/TDOT4
                  DTU1 = SQRT(A2**2 - 2.0*TD2/TDOT4) - A2
                  DTU1 = 1.01*MAX(DTU1,1.0D-10)
                  DTU = MIN(DTU1,DTU)
              END IF
          END IF
      END IF
*
*       Reset reference energy and generate new Stumpff coefficients.
      H0(IPAIR) = H(IPAIR)
   30 Z = -0.5D0*H(IPAIR)*DTU**2
      CALL STUMPF(IPAIR,Z)
      Z5 = SF(6,IPAIR)
      Z6 = SF(7,IPAIR)
      DT12 = ONE12*DTU*Z6
*
*       Convert to physical time units modified by Stumpff coefficients.
      STEP(I1) = (((((TDOT6*DT12 + TDOT5*Z5)*0.2*DTU + 0.5D0*TDOT4)*DTU
     &                     + TDOT3(IPAIR))*ONE6*DTU + TD2)*DTU + RI)*DTU
*
*       Ensure that regularized step is smaller than the c.m. step.
      IF (STEP(I1).GT.STEP(I).AND.HI.LT.0) THEN
          DTU = 0.5D0*DTU
          GO TO 30
      END IF
      DTAU(IPAIR) = DTU
*
*       See whether the KS slow-down procedure is activated.
      IMOD = KSLOW(IPAIR)
      IF (IMOD.GT.1) THEN
          ZMOD = FLOAT(ISLOW(IMOD))
          STEP(I1) = ZMOD*STEP(I1)
      END IF
*
*       Check diagnostics print option.
      IF (KZ(10).GE.3) THEN
          WRITE (6,40)  IPAIR, TIME+TOFF, H(IPAIR), RI, DTAU(IPAIR),
     &                  GI, STEP(I1), LIST(1,I1), IMOD
   40     FORMAT (3X,'KS MOTION',I6,2F10.4,1P,4E10.2,2I4)
      END IF
*
*       Employ special termination criterion in merger case.
      IF (NAME(I).LT.0) THEN
*       Terminate if apocentre perturbation > 0.25 (R > SEMI) or GI > 0.25.
          IF (HI.LT.0.0) THEN
*             SEMI = -0.5*BODY(I)/HI
*             ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
*             A0 = SEMI*(1.0 + SQRT(ECC2))/RI
*       Replace eccentricity calculation with typical value.
*             A0 = 1.5*SEMI/RI
*             GA = GI*A0*A0*A0
*             IF (GA.GT.0.25.AND.RI.GT.SEMI) IQ = .TRUE.
              IF (RI.GT.10*RMIN.AND.NNB0.GT.0.8*LIST(1,I)) IQ = .TRUE.
              IF (GI.GT.0.1.AND.RI.GT.RMIN) IQ = .TRUE.
              IF (GI.GT.0.01.AND.RI.GT.5.0*RMIN) IQ = .TRUE.
              IF (GI.GT.0.25) IQ = .TRUE.
*       Include extra condition for planet case.
              IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
                  IF (GI.GT.2.0D-04) IQ = .TRUE.
              END IF
          ELSE
              IF (TD2.GT.0.0.AND.(GI.GT.GMAX.OR.RI.GT.RMIN)) IQ = .TRUE.
          END IF
          IF (.NOT.IQ) GO TO 60
      END IF
*
*       Delay termination until end of block for large perturbations.
      IF (IQ) THEN
          DTR = TBLOCK - TIME
*         WRITE (6,45)  IPAIR, TTOT, GI, RI, DTR, STEP(I1)
*  45     FORMAT (' TERM TEST    KS T G R DTR DT  ',
*    &                           I4,F10.4,F7.3,1P,E10.2,2E9.1)
          IF (DTR.LT.STEP(I1)) GO TO 90
      END IF
*
*       Check standard termination criterion (suppress on IQ = .true.).
      IF (RI.GT.R0(IPAIR).AND.RI.GT.2.0*RMIN.AND..NOT.IQ) THEN
*       Include termination for rare tidal capture starting at pericentre.
          IF (KSTAR(I).LT.0.AND.RI.GT.5.0*RMIN) GO TO 90
*       Impose a limit using size of neighbour sphere (100*R > 0.80*RS).
          IF (RI.GT.8.0D-03*RS(I).AND.GI.GT.1.0D-04) GO TO 90
*       See whether termination can be delayed for sufficient perturbers.
          IF (NNB0.LT.0.80*LIST(1,I).AND.GI.LT.0.1) GO TO 60
*       Check updating of R0 for newly hardened binary orbit.
          IF (HI.LT.-ECLOSE) THEN
              SEMI = -0.5*BODY(I)/HI
              R0(IPAIR) = MAX(RMIN,2.0D0*SEMI) 
              R0(IPAIR) = MIN(R0(IPAIR),5.0*RMIN)
              GO TO 70
          END IF
*
*         IF (KZ(4).GT.0.AND.GI.GT.GPRINT(1)) THEN
*             DGAM = GI - G0
*             K = KZ(4)
*             DO 55 L = 2,K
*                 IF (GI.LT.GPRINT(L)) THEN
*                     CALL EVOLVE(IPAIR,L-1)
*                     GO TO 90
*                 END IF
*  55         CONTINUE
*             CALL EVOLVE(IPAIR,K)
*         END IF
          GO TO 90
      END IF
*
*       End integration cycle for hyperbolic motion.
   60 IF (HI.GE.0.0D0) THEN
          IF (RDOT*TD2.LT.0.0D0) THEN
*       Determine pericentre for hyperbolic two-body motion.
              SEMI = -0.5D0*BODY(I)/HI
              ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
              QPERI = SEMI*(1.0D0 - SQRT(ECC2))
              DMIN2 = MIN(DMIN2,QPERI)
*       MMadd: KZ(27) Tidal disruption diagnostics
              IF (KZ(19).GE.3.AND.KZ(27).EQ.3) THEN
                 IBHTD = 0
                 ITD = 0
                 ! Make sure the more massive body is allocated to BH
                 IF(BODY(I1).GT.BODY(I2).AND.BODY(I1).EQ.BODY1) THEN
                    IBHTD = I1
                    ITD = I2
                 ELSEIF(BODY(I2).GT.BODY(I1).AND.BODY(I2).EQ.BODY1) THEN
                    ITD = I1
                    IBHTD = I2
                 ENDIF
                 ! If we're including the BH
                 IF(IBHTD.GT.0) THEN
                    VINF = SQRT(2.0*HI)*VSTAR
                    KSBHTD = KSTAR(IBHTD)
                    KSTD = KSTAR(ITD)
                    ! Tidal radius
                    RTD = (BODY(IBHTD)/BODY(ITD))**0.3333 * RADIUS(ITD)
                    ! Factor x tidal radius for diagnostics
                    RTDF = 1.D3*RTD
                    ! tidal disruption
                    !WRITE(*,*) "About to check TD (hyper):",TIME,I1,I2
                    IF(QPERI.LE.RTD) THEN
                       ! write some diagnostics
!                       WRITE(*,*) "TIDAL DISRUPTION (hyper)", 
!     &                      NAME(IBHTD), NAME(ITD), QPERI,RTD
                       WRITE(57,*) TTOT, NAME(IBHTD),KSBHTD,BODY(IBHTD),
     &                      NAME(ITD),KSTD,BODY(ITD),RADIUS(ITD),
     &                      RTD,QPERI,SQRT(ECC2)
                       CALL FLUSH(57)
                       ! call CMBODY to do the encounter
                       CALL KSPERI(IPAIR)
                       KSPAIR = IPAIR
                       IQCOLL = -2
                       CALL CMBODY(QPERI,2)
                    ! tidal encounter without disruption
                    ELSEIF(QPERI.GT.RTD.AND.QPERI.LE.RTDF) THEN
!                       WRITE(*,*) "TIDAL ENCOUNTER (hyper)",
!     &                      NAME(IBHTD), NAME(ITD), QPERI,RTD
                       WRITE(58,*) TTOT, NAME(IBHTD),KSBHTD,BODY(IBHTD),
     &                      NAME(ITD),KSTD,BODY(ITD),RADIUS(ITD),
     &                      RTD,QPERI,SQRT(ECC2)
                       CALL FLUSH(58)
                    ENDIF
                 ENDIF ! if including BH
              ENDIF ! if extra KZ(27) TD ... end MMadd
*       Check optional tidal interaction or stellar collision.
*       MMadd: change type to include WDs
              IF (KZ(19).GE.3.AND.KSTAR(I).LE.14) THEN
                  VINF = SQRT(2.0*HI)*VSTAR
                  KS1 = KSTAR(I1)
                  KS2 = KSTAR(I2)
                  RX = MAX(RADIUS(I1),RADIUS(I2))
*       Determine maximum periastron factor for capture (VINF in km/sec).
                  IF (KZ(27).LE.2) THEN
                      RFAC = RPMAX2(RADIUS(I1),RADIUS(I2),BODY(I1),
     &                              BODY(I2),KS1,KS2,VINF)
                      RCAP = RFAC*RX
                  ELSE
                      DV = SQRT(2.0*HI)
*       Note that Quinlan & Shapiro function returns actual distance.
                      RCAP = RPMAX(BODY(I1),BODY(I2),VSTAR,DV,QPERI)
                  END IF
                  IF (QPERI.LT.5.0*RX) THEN
                      WRITE (54,54)  TTOT, NAME(I1), NAME(I2), KS1,
     &                               KS2, VINF, RCAP*SU, RX*SU, QPERI*SU
   54                 FORMAT (' CLOSE   T NAM K* VINF RCAP RX QP  ',
     &                                  F7.1,2I6,2I4,F6.2,3F6.1)
                  END IF
                  IF (QPERI.LT.RCAP) THEN
                      J1 = I1
                      IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                      FAC = 0.5*BODY(I)/BODY(J1)
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                      IF (KZ(27).LE.2) THEN
                          RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                      ELSE
                          RCOLL = 6.0*BODY(I)/CLIGHT**2
                      END IF
                      WRITE (55,58)  TTOT, IPAIR, NAME(I1), NAME(I2),
     &                               KS1, KS2, KSTAR(I), VINF
   58                 FORMAT (' RPMAX:    T KS NAM K* VINF ',
     &                                    F7.1,I5,2I6,3I4,F7.2)
                      WRITE (55,59)  SQRT(ECC2), HI, R(IPAIR), SEMI,
     &                               QPERI, BODY(I1), BODY(I2),
     &                               BODY(I)*ZMBAR
   59                 FORMAT (' E H R A QP BODY MT ',
     &                          F9.5,1P,6E10.2,0P,F6.1)
                      RI2 = 0.0
                      VI2 = 0.0
                      DO 61 K = 1,3
                          RI2 = RI2 + (X(K,I) - CMR(K))**2
                          VI2 = VI2 + XDOT(K,I)**2
   61                 CONTINUE
                      WRITE (55,62)  SQRT(RI2)/RC, SQRT(VI2)*VSTAR,
     &                               RHOD, RADIUS(I1)*SU, RADIUS(I2)*SU,
     &                               RCAP, RADIUS(J1)/QPERI, RCOLL/QPERI
   62                 FORMAT (' r/RC V* <C> R* RCAP R1/QP RCOLL/QP ',
     &                          2F5.1,3F6.1,3F5.1)
                      CALL FLUSH(55)
                      IF (QPERI.LT.RCOLL) THEN
*       Obtain KS variables at pericentre before merging into one body.
                          CALL KSPERI(IPAIR)
                          KSPAIR = IPAIR
                          IQCOLL = -2
      WRITE (6,600)  NSTEPU, QPERI
  600 FORMAT (' CALL CMBODY   # QP ',I11,1P,E10.2)
                          CALL CMBODY(QPERI,2)
                      ELSE IF (KSTAR(I).GE.0.AND.KZ(27).GT.0) THEN
                          CALL KSTIDE(IPAIR,QPERI)
                      END IF
                  END IF
*       Check options for artificial collisions.
              ELSE IF (KZ(27).EQ.-1.AND.KZ(13).LT.0) THEN
                  RFAC = 2.0
                  IF (QPERI.LT.RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                      J1 = I1
                      IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                      FAC = 0.5*BODY(I)/BODY(J1)
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                      RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                      IF (QPERI.LT.RCOLL) THEN
                          CALL TOUCH(IPAIR,I1,I2,RCOLL)
                      END IF
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Check perturbation threshold (H < 0 & GAMMA > GMAX).
*     IF (KZ(4).EQ.0.OR.G0.LT.GMAX) GO TO 70
*
*     K = KZ(4)
*     DO 65 L = 1,K
*         IF ((G0 - GPRINT(L))*(GI - GPRINT(L)).LE.0.0) THEN
*             IF (L.EQ.1) THEN
*                 W1 = -0.5*BODY(I)/HI
*                 W2 = W1*BODYIN
*                 TK = TWOPI*W1*SQRT(W2)
*             END IF
*
*       Estimate smallest permitted output interval at new level.
*             DTCRIT = TK*ORBITS(L)
*             IF (TIME - TLASTB(L).LT.DTCRIT) GO TO 70
*             DGAM = GI - G0
*             CALL EVOLVE(IPAIR,L)
*             GO TO 70
*         END IF
*  65 CONTINUE
*
*       Check for partial reflection during approach (NB! only IMOD = 1).
*  70 IF (GI.LT.GMIN.AND.TD2.LT.0.0D0) THEN
*       Skip apocentre position itself.
*         IF (RDOT.LT.0.0D0.AND.IMOD.EQ.1) THEN
*             IF (KZ(25).GT.0) CALL FREEZE(IPAIR)
*             GO TO 100
*         END IF
*     END IF
*
*       Determine new perturbers for binary at apocentre turning point.
   70 IF (RDOT*TD2.GE.0.0D0) GO TO 100
      SEMI = -0.5D0*BODY(I)/HI
*
*       Check minimum two-body separation just after pericentre.
      IF (RDOT.LT.0.0D0) THEN
*       Obtain pericentre by Mikkola's algorithm (GAMMA < 0.001).
          IF (GI.LT.0.001) THEN
              CALL PERI(UI,UIDOT,RI,BODY(I1),BODY(I2),QPERI)
          ELSE
              QPERI = RI
          END IF
          DMIN2 = MIN(DMIN2,QPERI)
*     MMadd: KZ(27) Tidal disruption diagnostics
          IF (KZ(19).GE.3.AND.KZ(27).EQ.3) THEN
             IBHTD = 0
             ITD = 0
             ! Make sure the more massive body is allocated to BH
             IF(BODY(I1).GT.BODY(I2).AND.BODY(I1).EQ.BODY1) THEN
                IBHTD = I1
                ITD = I2
             ELSEIF(BODY(I2).GT.BODY(I1).AND.BODY(I2).EQ.BODY1) THEN
                ITD = I1
                IBHTD = I2
             ENDIF
             ! If we're including the BH
             IF(IBHTD.GT.0) THEN
                VINF = SQRT(2.0*HI)*VSTAR
                ECC2 = (1.0 - RI/SEMI)**2 +
     &               TDOT2(IPAIR)**2/(BODY(I)*SEMI)
                KSBHTD = KSTAR(IBHTD)
                KSTD = KSTAR(ITD)
                ! Tidal radius
                RTD = (BODY(IBHTD)/BODY(ITD))**0.3333 * RADIUS(ITD)
                ! Factor x tidal radius for diagnostics
                RTDF = 1.D3*RTD
                ! tidal disruption
                !WRITE(*,*) "About to check TD (ecc):",TIME,I1,I2
                IF(QPERI.LE.RTD) THEN
                   ! write some diagnostics
 !                  WRITE(*,*) "TIDAL DISRUPTION (ecc)", 
 !    &                  NAME(IBHTD), NAME(ITD), QPERI
                   WRITE(57,*) TTOT, NAME(IBHTD),KSBHTD,BODY(IBHTD),
     &                  NAME(ITD),KSTD,BODY(ITD),RADIUS(ITD),
     &                  RTD,QPERI,SQRT(ECC2) 
                   CALL FLUSH(57)
                   ! call CMBODY to do the encounter
                   CALL KSPERI(IPAIR)
                   KSPAIR = IPAIR
                   IQCOLL = -2
                   CALL CMBODY(QPERI,2)
                   ! tidal encounter without disruption
                ELSEIF(QPERI.GT.RTD.AND.QPERI.LE.RTDF) THEN
!                   WRITE(*,*) "TIDAL ENCOUNTER (ecc)",
!     &                  NAME(IBHTD), NAME(ITD), QPERI
                   WRITE(58,*) TTOT,NAME(IBHTD),KSBHTD,BODY(IBHTD),
     &                  NAME(ITD),KSTD,BODY(ITD),RADIUS(ITD),
     &                  RTD,QPERI,SQRT(ECC2)
                   CALL FLUSH(58)
                ENDIF
                ! check the inspiral case if compact objects
                IF(KSTD.GE.10) THEN
                   CALL TIDES3(QPERI,BODY(I1),BODY(I2),
     &                  VSTAR,H(IPAIR),SQRT(ECC2),DE)
                   ! DH = DE/MU
                   DH = DE*(BODY(I1)+BODY(I2))/(BODY(I1)*BODY(I2))
                   TINSP = ABS(H(IPAIR)/DH)*SQRT(SEMI**3/BODY(IPAIR))
!                   WRITE(*,*), "COMPACT BIN, Tinsp = ",
!     &                  NAME(I1),NAME(I2),SEMI,TINSP
                   IF(QPERI.LT.1.d-5.AND.(TT0T-IDINT(TTOT)).LT.0.01)THEN
                      WRITE (*,*) "COMPACT BINARY",TTOT,NAME(IBHTD),
     &                             NAME(ITD),TINSP
                   ENDIF
                   IF(TINSP.LT.1.D0) THEN
                      ! call CMBODY to do the encounter                                                                                                                                   
                      WRITE(*,*) "COMPACT INSPIRAL:",TTOT,
     &                     NAME(IBHTD),KSBHTD,BODY(IBHTD),
     &                     NAME(ITD),KSTD,BODY(ITD),RADIUS(ITD),
     &                     QPERI,SEMI, SQRT(ECC2)
                      CALL KSPERI(IPAIR)
                      KSPAIR = IPAIR
                      IQCOLL = -2
                      CALL CMBODY(QPERI,2)
                      
                   ENDIF !inspiral
                ENDIF! check inspiral
             ENDIF              ! if including BH
          ENDIF                 ! if extra KZ(27) TD ... end MMadd
*
*       Check optional tidal interaction or stellar collision (skip merger).
          IF (KZ(19).GE.3.AND.KSTAR(I).LE.14.AND.NAME(I).GT.0) THEN
              RFAC = 5.0
              IF (KZ(27).LE.2) THEN
                  IF (KZ(27).EQ.1) RFAC = 4.0
                  RX = RFAC*MAX(RADIUS(I1),RADIUS(I2))
              ELSE
                  RX = RPMIN(BODY(I1),BODY(I2),VSTAR,HI,QPERI)
              END IF
              IF (QPERI.LT.RX) THEN
                  J1 = I1
                  IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                  FAC = 0.5*BODY(I)/BODY(J1)
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                  IF (KZ(27).LE.2) THEN
                      RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                  ELSE
                      RCOLL = 6.0*BODY(I)/CLIGHT**2
                  END IF
                  IF (QPERI.LT.RCOLL) THEN
*       Obtain KS variables at pericentre before merging into one body.
                     WRITE(*,*) "QPERI<RCOLL (ecc)", QPERI, RCOLL
                     CALL KSPERI(IPAIR)
                     KSPAIR = IPAIR
                     IQCOLL = -2
                     CALL CMBODY(QPERI,2)
                  ELSE IF (KSTAR(I).GE.0) THEN
*       Distinguish between sequential, standard and GR circularization.
                      IF (KZ(27).EQ.1) THEN
                          ICIRC = 1
                          TC = 0.0
                      ELSE IF (KZ(27).EQ.2.AND.KSTAR(I).LT.14) THEN
                          ECC2 = (1.0 - RI/SEMI)**2 +
     &                                    TDOT2(IPAIR)**2/(BODY(I)*SEMI)
                          ECC = SQRT(ECC2)
                          ICIRC = 0
                          CALL TCIRC(QPERI,ECC,I1,I2,ICIRC,TC)
                      ELSE
                          ICIRC = 1
                          TC = 0.0
                      END IF
!MMcommenting this                      IF (KSTAR(I).GE.10) ICIRC = 0
*       Skip tidal effects for circularization time above 100 Myr (07/08).
                      IF (ICIRC.GT.0.AND.KZ(27).GT.0.AND.
     &                    TC.LT.100.0) THEN
                         !WRITE(*,*) "QPERI>RCOLL (ecc)", QPERI, RCOLL
                         CALL KSTIDE(IPAIR,QPERI)
                      END IF
                  END IF
              END IF
*       Check for perturbed spiral or chaos case (skip collision).
              IF (KSTAR(I).EQ.-2.AND.IPHASE.EQ.0) THEN
                  CALL SPIRAL(IPAIR)
              ELSE IF (KSTAR(I).EQ.-1.AND.IPHASE.EQ.0) THEN
                  CALL KSTIDE(IPAIR,QPERI)
              END IF
*       Check options for artificial collisions.
          ELSE IF (KZ(27).EQ.-1.AND.KZ(13).LT.0) THEN
              RFAC = 2.0
              IF (QPERI.LT.RFAC*MAX(RADIUS(I1),RADIUS(I2))) THEN
                  J1 = I1
                  IF (RADIUS(I2).GT.RADIUS(I1)) J1 = I2
                  FAC = 0.5*BODY(I)/BODY(J1)
*       Adopt collision criterion of Kochanek (Ap.J. 385, 604, 1992).
                  RCOLL = 1.7*FAC**0.3333*RADIUS(J1)
                  IF (QPERI.LT.RCOLL) THEN
                      CALL TOUCH(IPAIR,I1,I2,RCOLL)
                  END IF
              END IF
          END IF
          GO TO 100
      END IF
*
*       Save maximum separation of persistent binary.
      RMAX = MAX(RMAX,RI)
*
*       Check binary reference radius or merger stability criterion.
      IF (NAME(I).GT.0) THEN
*       Update termination length scale in case of initial soft binary.
          EB = BODY(I1)*BODY(I2)*HI*BODYIN
          IF (EB.LT.EBH) R0(IPAIR) = MAX(RMIN,2.0*SEMI)
      ELSE 
          ECC2 = (1.0 - RI/SEMI)**2 + TDOT2(IPAIR)**2/(BODY(I)*SEMI)
          ECC = SQRT(ECC2)
          RP = SEMI*(1.0 - ECC)*(1.0 - 2.0*GI)
*       Find merger index.
          IM = 0
          DO 72 K = 1,NMERGE
              IF (NAMEM(K).EQ.NAME(I)) IM = K
   72     CONTINUE
*       Exclude inner planets from the general stability test.
          IF (MIN(BODY(I1),BODY(I2)).LT.0.05*BODYM) THEN
              IF (RP.LT.R0(IPAIR)) GO TO 90
          END IF
*       Include optional Kozai diagnostics (large EMAX) & circularization.
          IF (KZ(42).GT.0) THEN
              CALL KOZAI(IPAIR,IM,ECC,SEMI,ITERM)
*       Check merger termination and perform circularization or collision.
              IF (ITERM.GT.0) THEN
                  IPHASE = 7
                  KSPAIR = IPAIR
                  CALL RESET
                  IF (ITERM.GT.1) CALL CIRC
                  GO TO 100
              END IF
          END IF
*       Assess the stability inside critical pericentre (safety factor 1.04).
          IF (RP.LT.1.04*R0(IPAIR)) THEN
*       Note: assessment needs to use same eccentricity as for acceptance.
              CALL ASSESS(IPAIR,IM,ECC,SEMI,ITERM)
              IF (ITERM.GT.0) THEN
                  INSTAB = INSTAB + 1
                  GO TO 90
              END IF
          END IF
*       Check possible eccentricity modulation or t_circ update.
          IF (IM.GT.0.AND.(TIME.GT.TMDIS(IM).OR.
     &        TMDIS(IM).GT.1.0D+06)) THEN
              IF (KZ(27).EQ.2) THEN
                  CALL ECCMOD(I,ITERM)
                  IF (ITERM.GT.0) THEN
*                     WRITE (6,76)  RP, R0(IPAIR)
*  76                 FORMAT (' ECCMOD TERM    RP R0 ',1P,2E10.2)
                      GO TO 90
                  END IF
*       Consider both inner and possible outer circularizing binary.
                  DO 78 K = 1,NCHAOS
                      IF (NAMEC(K).EQ.NZERO - NAMEM(IM).AND.
     &                    KSTARM(IM).EQ.-2) THEN
*       Update unperturbed binary if T - TOSC > 10 Myr (cf. IMPACT & DECIDE).
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) GO TO 90
                      END IF
                      IF (NAMEC(K).EQ.NAMEG(IM).AND.
     &                    KSTARM(IM).EQ.-2) THEN
                          IF ((TIME - TOSC(K))*TSTAR.GT.10.0) GO TO 90
                      END IF
*       Note: perturbed binary is treated if pericentre before next IMPACT.
   78             CONTINUE
              END IF
          END IF
      END IF
*
*       Produce diagnostics for any circularizing perturbed binary.
      IF (KSTAR(I).EQ.-2.AND.GI.GT.0.01) THEN
          ECC = RI/SEMI - 1.0
          QPS = SEMI*(1.0 - ECC)/MAX(RADIUS(I1),RADIUS(I2))
          ZM = BODY(I1)/BODY(I2)
          WRITE (21,80)  TTOT, NAME(I1), NAME(I2), KSTAR(I1), KSTAR(I2),
     &                   LIST(1,I1), QPS, ZM, GI, ECC, SEMI
   80     FORMAT (' PERT SPIRAL    T NAM K* NP QP/S M1/M2 G E A ',
     &                             F11.4,2I6,3I3,2F5.1,2F8.4,1P,E10.2)
          CALL FLUSH(21)
      END IF
*
*       See whether KS slow-down procedure should be (re)-checked (no Chaos).
      IF (KZ(26).GT.0.AND.KSTAR(I).GE.0) THEN
          KMOD = RANGE*GMIN/MAX(GI,1.0D-10)
          IF (KMOD.GT.1.OR.IMOD.GT.1) THEN
              CALL KSMOD(IPAIR,KMOD)
              IF (KMOD.LT.0) GO TO 100
              GO TO 82
          END IF
      END IF
*
*       Set approximate value of next period with perturbation included.
      TK = TWOPI*SEMI*SQRT(SEMI*BODYIN)*(1.0 + GI)
      IF (IMOD.GT.1) THEN
          TK = ZMOD*TK
      END IF
*
*       Use old perturber list if next apocentre is before the c.m. step.
      IF (TIME + TK.LT.T0(I) + STEP(I)) THEN
          GO TO 100
      END IF
*
*       Select new perturbers (J = N adopted for unperturbed Chaos).
   82 CALL KSLIST(IPAIR)
*
*       Check rectification of chaotic spiral at start of unperturbed motion.
      IF (KSTAR(I).EQ.-2.AND.LIST(1,I1).EQ.0) THEN
          DMR = 0.D0
          CALL CHRECT(IPAIR,DMR)
          IF (IPHASE.LT.0) GO TO 100
      ELSE
          CALL KSRECT(IPAIR)
      END IF
*
*       See whether a massive BH subsystem can be selected.
      IF (KZ(11).NE.0.AND.NCH.EQ.0.AND.BODY(I).GT.10.0*BODYM.AND.
     &   SEMI.LT.7.0*RMIN.AND.LIST(1,I1).LE.10.AND.NAME(I).GT.0) THEN
*    &   SEMI.LT.RMIN.AND.LIST(1,I1).LE.5.AND.NAME(I).GT.0) THEN
*
*     IF (SEMI.GT.0.1*RMIN) GO TO 88
*       Check optional BH condition (prevents mass-loss complications).
          IF (KZ(11).LE.-2) THEN
              IF (KSTAR(I1).NE.14.OR.KSTAR(I2).NE.14) GO TO 88
          END IF
*
          RCR2 = RMIN2*BODY(I)/(10.0*BODYM)
          NNB1 = LIST(1,I1) + 1
          RX2 = 1000.0
          JCLOSE = 0
          DO 85 L = 2,NNB1
              J = LIST(L,I1)
              RIJ2 = (X(1,I) - X(1,J))**2 + (X(2,I) - X(2,J))**2 +
     &                                      (X(3,I) - X(3,J))**2
              RD = (X(1,I)-X(1,J))*(XDOT(1,I)-XDOT(1,J)) +
     &             (X(2,I)-X(2,J))*(XDOT(2,I)-XDOT(2,J)) +
     &             (X(3,I)-X(3,J))*(XDOT(3,I)-XDOT(3,J))
              IF (RIJ2.LT.RCR2.AND.RD.LT.0.0.OR.
     &           (RIJ2.LT.4.0*SEMI**2)) THEN
                  IF (RIJ2.LT.RX2) THEN
                      RX2 = RIJ2
                      RRD = RD
                      JCLOSE = J
                  END IF
              END IF
   85     CONTINUE
          IF (JCLOSE.GT.0) THEN
              IF (NAME(JCLOSE).LE.0) GO TO 88
              RX = SQRT(RX2)
*       Limit energy of triple system (< 50*EBH) using radial velocity.
              ZMU = BODY(I)*BODY(JCLOSE)/(BODY(I) + BODY(JCLOSE))
              EBT = EB + ZMU*(0.5*(RRD/RX)**2-(BODY(I)+BODY(JCLOSE)/RX))
              IF (EBT.GT.50.0*EBH) GO TO 88
              WRITE (6,86)  TTOT, NAME(JCLOSE), LIST(1,I1), STEP(I),
     &                      STEP(JCLOSE), SEMI, RX, GI, EBT
   86         FORMAT (' NEW CHAIN   T NMJ NP STEPI STEPJ A RIJ GI EBT ',
     &                              F9.3,I6,I4,1P,6E10.2)
              CALL FLUSH(6)
*       Initiate chain regularization directly (B-B or B-S: see IMPACT).
              JCOMP = JCLOSE
              JCMAX = 0
              KSPAIR = IPAIR
              IPHASE = 8
              EBCH0 = EB
*       Distinguish between case of single and binary intruder.
              IF (JCLOSE.LE.N) THEN
                  KS2 = 0
              ELSE
                  KS2 = JCLOSE - N
*       Reduce c.m. body location and pair index if IPAIR terminated first.
                  IF (KS2.GT.IPAIR) JCLOSE = JCLOSE - 1
                  IF (KS2.GT.IPAIR) KS2 = KS2 - 1
              END IF
*       Initialize new ARchain.
              CALL DELAY(1,KS2)
              GO TO 100
          END IF
      END IF
*
*       Check optional search criterion for multiple encounter or merger.
   88 IF (KZ(15).GT.0.AND.STEP(I).LT.DTMIN) THEN
          CALL IMPACT(I)
      END IF
!     MMadd Morgan adding own impact routine
      IF(KZ(45) .GE. 3 ) CALL IMPMM(I)
      GO TO 100
*
*       Terminate regularization of current pair (IPAIR set in KSPAIR).
   90 KSPAIR = IPAIR
*       Set indicator for calling KSTERM in MAIN (permits phase overlay).
      IPHASE = 2
*       Check case of hierarchical binary.
      IF (NAME(I).LT.0) IPHASE = 7
*
  100 RETURN
*
      END
