      SUBROUTINE EXPEL(J1,J2,ICASE)
*
*       Common envelope stage of interacting stars & BH.
*       ------------------------------------------------
*
      INCLUDE 'common6.h'
      COMMON/MODES/  EB0(NTMAX),ZJ0(NTMAX),ECRIT(NTMAX),AR(NTMAX),
     &               BR(NTMAX),EOSC(4,NTMAX),EDEC(NTMAX),TOSC(NTMAX),
     &               RP(NTMAX),ES(NTMAX),CM(2,NTMAX),IOSC(NTMAX),
     &               NAMEC(NTMAX)
      REAL*8 M01,M1,R1,AJ1,M02,M2,R2,AJ2,SEP,MI(2)
      REAL*8 TSCLS(20),LUMS(10),GB(10),TM,TN,LUM1,LUM2,MC1,MC2,RCC
      REAL*8 JSPIN1,JSPIN2,MENV,RENV,K2
      CHARACTER*8  WHICH1
      LOGICAL COALS
*
*
*       Define global indices such that body #I1 is giant-like.
      IF(KSTAR(J1).GE.2.AND.KSTAR(J1).LE.9.AND.KSTAR(J1).NE.7)THEN
          I1 = J1
          I2 = J2
      ELSE
          I1 = J2
          I2 = J1
      ENDIF
*
*       Switch to #I2 if #I1 is a BH (standard case not affected).
      IF (KSTAR(I1).EQ.14) THEN
          II = I1
          I1 = I2
          I2 = II
      END IF
      !MMadd, if they're both BHs, still want I2 to be more massive
      IF(KSTAR(I1).EQ.14.AND.KSTAR(I2).EQ.14) THEN
         IF(BODY(I1).GT.BODY(I2))THEN
            II=I1
            I1=I2
            I2=II
         ENDIF
      ENDIF

*
*       Update the stars to latest previous time.
      TEV1 = MAX(TEV0(I1),TEV0(I2))
*
*       Specify basic parameters for both stars.
      M01 = BODY0(I1)*ZMBAR
      M1 = BODY(I1)*ZMBAR
      MC1 = 0.D0
      AJ1 = TEV1*TSTAR - EPOCH(I1)
      JSPIN1 = SPIN(I1)*SPNFAC
      KW1 = KSTAR(I1)
      M02 = BODY0(I2)*ZMBAR
      M2 = BODY(I2)*ZMBAR
      MC2 = 0.D0
      AJ2 = TEV1*TSTAR - EPOCH(I2)
      JSPIN2 = SPIN(I2)*SPNFAC
      KW2 = KSTAR(I2)
      ITER = 0
*
*       Save current semi-major axis.
      IPAIR = KSPAIR
      I = N + IPAIR
      SEMI0 = -0.5d0*BODY(I)/H(IPAIR)
      SEP = SEMI0*SU
      ECC2 = (1.D0-R(IPAIR)/SEMI0)**2+TDOT2(IPAIR)**2/(SEMI0*BODY(I))
      ECC = SQRT(ECC2)
      ECC0 = ECC
      DM2 = 0.0
      EP2 = EPOCH(I2)
      AGE2 = AJ2
      RP1 = R(IPAIR)
      ISKIP = 0
*
*       Include special disruption treatment for BH + star (only #11 < 0).
      IF (MAX(KSTAR(I1),KSTAR(I2)).EQ.14) THEN
          !DM1 = 0.5*BODY(I1)*ZMBAR
          !MMadd: make dm = mstar
          DM1 = BODY(I1)*ZMBAR
        
*       Conserve total mass and specific energy (H = const here).
          M1 = M1 - DM1
          M2 = M2 + DM1
          R1 = RADIUS(I1)*SU
          R2 = RADIUS(I2)*SU
*       Transform to positive radial velocity and obtain global coordinates.
          ITER = 0
    2     JPAIR = 0
          CALL KSAPO(JPAIR)   ! note JPAIR = KSPAIR on return.
*       Note zero argument which advances eccentric anomaly by 0.02.
          T0(I1) = TBLOCK
          ITER = ITER + 1
          IF (ITER.GT.10) WRITE (6,70)  ITER, R(IPAIR), TDOT2(IPAIR)
   70     FORMAT (' KSAPO!    IT R TD2 ',I4,1P,2E10.2)
          IF (ITER.LT.50.AND.R(IPAIR).LT.RMIN) GO TO 2
          CALL RESOLV(IPAIR,1)
          ! MMadd: set COALS = true
          COALS = .TRUE.
          !COALS = .FALSE.
          !IF(KSTAR(I2).EQ.14.AND.KSTAR(I1).GE.13) COALS = .TRUE.
          SEMI = SEMI0
          ISKIP = 1
          !KW1 = 15 ! send disrupted object as massless ghost
*       Skip common envelope part for BH + star interaction
          GO TO 4
      END IF
*
*       Perform common envelope evolution (note: SEP in SU).
    1 CALL COMENV(M01,M1,MC1,AJ1,JSPIN1,KW1,
     &            M02,M2,MC2,AJ2,JSPIN2,KW2,ECC,SEP,COALS)
*
*       Obtain consistent radii for the stars (skip #I2 on coalescence).
      CALL star(KW1,M01,M1,TM,TN,TSCLS,LUMS,GB,ZPARS)
      CALL hrdiag(M01,AJ1,M1,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &            R1,LUM1,KW1,MC1,RCC,MENV,RENV,K2)
      IF(KW1.NE.KSTAR(I1))THEN
         WRITE(38,*)' EXPEL TYPE CHANGE1 ',KSTAR(I1),KW1
      ENDIF
      IF(COALS)THEN
          KW2 = KSTAR(I2)
          R2 = 0.d0
          LUM2 = 0.d0
          JSPIN2 = 0.d0
      ELSE
          CALL star(KW2,M02,M2,TM,TN,TSCLS,LUMS,GB,ZPARS)
          CALL hrdiag(M02,AJ2,M2,TM,TN,TSCLS,LUMS,GB,ZPARS,
     &                R2,LUM2,KW2,MC2,RCC,MENV,RENV,K2)
          IF(KW2.NE.KSTAR(I2))THEN
             WRITE(38,*)' EXPEL TYPE CHANGE2 ',KSTAR(I2),KW2
          ENDIF
      END IF
*
*       Check coalescence condition again (bug 7/1/09).
      IF (SEP.LT.R1.OR.SEP.LT.R2) THEN
          COALS = .TRUE.
      END IF
*
*       Skip small mass loss unless coalescence (May 2003).
      DMS = M1 + M2 - BODY(I)*SMU
      IF (ABS(DMS).LT.1.0D-10.AND..NOT.COALS) THEN
*         WRITE (77,3)  TIME, ECC, DMS, SEMI0*SU-SEP
*   3     FORMAT (' FAIL!    T ECC DMS DSEP ',F12.4,F8.4,1P,2E10.2)
*         CALL FLUSH(77)
          IPHASE = 0
          GO TO 60
      END IF
*
*       Copy new values of radius, luminosity & semi-major axis.
      RADIUS(I1) = R1/SU
      RADIUS(I2) = R2/SU
      ZLMSTY(I1) = LUM1
      ZLMSTY(I2) = LUM2
      SPIN(I1) = JSPIN1/SPNFAC
      SPIN(I2) = JSPIN2/SPNFAC
*       Accept negative semi-major axis for fly-by (bug fix 12/08).
      SEMI = SEP/SU
*
*       Set reduced mass & binding energy and update masses & time T0(J1).
    4 ZMU0 = BODY(I1)*BODY(I2)/BODY(I)
      HI = H(IPAIR)
      BODY0(I1) = M01/ZMBAR
      BODY0(I2) = M02/ZMBAR
      T0(J1) = TIME
*       Add current energy (final state: coalescence, binary or disruption).
      ECOLL = ECOLL + ZMU0*HI
      EGRAV = EGRAV + ZMU0*HI
*
*       Distinguish between coalescence and surviving binary.
      EPOCH(I1) = TEV1*TSTAR - AJ1
      IF(COALS)THEN
          MI(1) = M1
          MI(2) = M2
          BODY0(I2) = 0.D0
*       Check for TZ object formed by CE evolution.
          IF(KSTAR(I2).GE.13.AND.KW1.GE.13)THEN
              NTZ = NTZ + 1
              WRITE (6,5)  KW1, M1, M2
    5         FORMAT (' NEW TZ    K* M1 M2 ',I4,2F7.2)
          ENDIF
          CALL COAL(IPAIR,KW1,KW2,MI)
      ELSE
*
*       Update evolution times and TMDOT.
          EPOCH(I2) = TEV1*TSTAR - AJ2
          TEV(I1) = TEV1
          TEV0(I1) = TEV(I1)
          TEV(I2) = TEV1
          TEV0(I2) = TEV(I2)
          TMDOT = MIN(TMDOT,TEV1)
*       Shrink the orbit in case of no coalescence.
          BODY(I1) = M1/ZMBAR
          BODY(I2) = M2/ZMBAR
*       Form mass loss & new binary mass and re-define binding energy.
          DM = BODY(I) - (BODY(I1) + BODY(I2))
*       Save mass loss for potential energy correction after ENFORCED CE.
          DM2 = DM2 + DM
          ZMASS = ZMASS - DM
          BODY(I) = BODY(I) - DM
          H(IPAIR) = -0.5d0*BODY(I)/SEMI
*
*       Subtract final binding energy of surviving binary.
          ZMU = BODY(I1)*BODY(I2)/BODY(I)
          ECOLL = ECOLL - ZMU*H(IPAIR)
*         EGRAV = EGRAV - ZMU*H(IPAIR)
*
          IF (ISKIP.EQ.0) THEN
*       Set argument for new pericentre velocity (need to have H*A*(1-E)).
              RX = R(IPAIR)/(1.0 - ECC)
*       Modify KS variables consistent with new eccentricity and energy H.
              CALL EXPAND(IPAIR,RX)
              CALL RESOLV(IPAIR,1)
          ELSE
*       Fine-tune the variables (especially for hyperbolic pericentre).
              CALL RESOLV(IPAIR,1)
          END IF
*
          IF (ECC0.LT.1.0) THEN
              WHICH1 = ' BINARY '
          ELSE
              WHICH1 = ' HYPERB '
              NHYPC = NHYPC + 1
          END IF
          IF (ISKIP.GT.0)  THEN
              WHICH1 = 'DISRUPT '
              NDISR = NDISR + 1
*       Restore mass of accreted star for diagnostics (from 23/5/14).
              M1 = M1+DM1
              WRITE (24,*)  TIME+TOFF, NDISR, NAME(I1), KSTAR(I1),
     &                      ECC, MIN(M1,M2), MAX(M1,M2), SEMI
    6         FORMAT (' DISRUPT1    T NDISR NM K* E M1 M2 SEMI ',
     &                              F8.1,I5,I7,I4,F10.6,2F6.1,1P,E10.2)
              CALL FLUSH(24)
          END IF
*
          WRITE (6,10)  WHICH1, NAME(I1), NAME(I2), KSTAR(I1),KSTAR(I2),
     &                  KW1, KW2, M1, M2, DM*ZMBAR, ECC0, ECC, R1, R2,
     &                  RP1, SEMI0*SU, SEMI*SU
   10     FORMAT (A8,'CE    NAM K0* K* M1 M2 DM E0 E R1 R2 RP A0 A ',
     &               2I6,4I3,ES15.3,ES15.3,ES15.3,2F10.6,2F7.1,1P3E9.1)
*
*       Check common envelope condition again after circularization (09/08).
          IF (ECC0.GT.0.001D0.AND.ECC.LE.0.001D0) THEN
              Q1 = M1/M2
              Q2 = 1.D0/Q1
              RL1 = RL(Q1)
              RL2 = RL(Q2)
              IF (R1.GT.RL1*SEP.OR.R2.GT.RL2*SEP) THEN
                  ITER = ITER + 1
                  IF (ITER.EQ.1) THEN
                      WRITE (6,8)  RL1*SEP, RL2*SEP
    8                 FORMAT (' ENFORCED CE    RL*SEP ',1P,2E10.2)
                      ECC0 = ECC
                      SEMI0 = SEMI
                      SEP = SEMI*SU
*       Correct for mass loss in case of coalescence which ignores DM2.
                      IF (DM*SMU.LT.0.05) THEN
                          CALL FICORR(I,DM)
                      ELSE
                          CALL FCORR(I,DM,0)
                      END IF
                      DM2 = 0.0
                      GO TO 1
                  END IF
              END IF
          END IF
*
*       Copy mass loss (one or two COMENV) and new types.
          IF (ITER.GT.0) DM = DM2
          KSTAR(I1) = KW1
          KSTAR(I2) = KW2
*
*       Save perturber list for force corrections and polynomials.
          NNB = LIST(1,I1)
          DO 11 L = 1,NNB+1
              ILIST(L) = LIST(L,I1)
   11     CONTINUE
*
          write(*,*) "testing for massless ghost"
*       Ensure rare case of massless remnant will escape at next output.
          IF (KW1.EQ.15.OR.KW2.EQ.15) THEN
            CALL KSTERM
            IPHASE = -1
            I = IFIRST
            IF (BODY(I).GT.0.0D0) I = IFIRST + 1
            T0(I) = TADJ + DTADJ
            STEP(I) = 1.0D+06
            RI = SQRT(X(1,I)**2 + X(2,I)**2 + X(3,I)**2)
            VI = SQRT(XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2)
            DO 12 K = 1,3
               X0(K,I) = MIN(1.0d+04+X(K,I),1000.0*RSCALE*X(K,I)/RI)
               X(K,I) = X0(K,I)
               X0DOT(K,I) = SQRT(0.004*ZMASS/RSCALE)*XDOT(K,I)/VI
               XDOT(K,I) = X0DOT(K,I)
               F(K,I) = 0.0D0
               FDOT(K,I) = 0.0D0
               D2(K,I) = 0.0D0
               D3(K,I) = 0.0D0
   12       CONTINUE
            WRITE (6,13)  NAME(I), KSTAR(I), BODY(I)*ZMBAR
   13       FORMAT (' MASSLESS GHOST    NAM K* M ',I6,I4,1P,E10.2)
*
*       Include special treatment for velocity kick of KS binary.
          ELSE IF (KW1.GE.10.OR.ISKIP.GT.0) THEN
*
              IF (ECC.LE.0.001D0) KSTAR(N+IPAIR) = 10
*       Re-initialize the KS solution with new perturber list.
              CALL KSLIST(IPAIR)
              CALL KSPOLY(IPAIR,1)
*       Evaluate binary disruption velocity and introduce velocity kick.
              VI2 = XDOT(1,I)**2 + XDOT(2,I)**2 + XDOT(3,I)**2
              KSTAR(I1) = -KSTAR(I1)
              IF (ISKIP.GT.0) THEN
*       Ensure binary velocity kick by setting BH index (not saved).
                  KW1 = 14
              END IF
              CALL KICK(IPAIR,0,KW1,DM)
*
*       Include potential and kinetic energy corrections due to mass loss.
              CALL POTI(I,POTJ)
              ECDOT = ECDOT + DM*POTJ + 0.5d0*DM*VI2
*
*       See whether tidal terms should be included (standard or scaled).
              IF (KZ(14).GT.0) THEN
                  IF (KZ(14).LE.2) THEN
                      ECDOT = ECDOT - 0.5d0*DM*(TIDAL(1)*X(1,I)**2 +
     &                                          TIDAL(3)*X(3,I)**2)
                  ELSE
                      BODY(I) = BODY(I) + DM
                      CALL XTRNLV(I,I)
                      ECDOT = ECDOT + HT
                      BODY(I) = BODY(I) - DM
                      CALL XTRNLV(I,I)
                      ECDOT = ECDOT - HT
                  END IF
              END IF
*
*       Form new c.m. variables.
              T0(I) = TIME
              DO 20 K = 1,3
                  X(K,I) = (BODY(I1)*X(K,I1) + BODY(I2)*X(K,I2))/BODY(I)
                  XDOT(K,I) = (BODY(I1)*XDOT(K,I1) +
     &                         BODY(I2)*XDOT(K,I2))/BODY(I)
                  X0(K,I) = X(K,I)
                  X0DOT(K,I) = XDOT(K,I)
   20         CONTINUE
*
*       Obtain new F & FDOT and time-steps for neighbours.
              DO 30 L = 2,NNB+1
                  J = ILIST(L)
                  DO 25 K = 1,3
                      X0DOT(K,J) = XDOT(K,J)
   25             CONTINUE
                  CALL FPOLY1(J,J,0)
                  CALL FPOLY2(J,J,0)
   30         CONTINUE
              TPREV = TIME - STEPX
*
*       Include safety condition to avoid ejection of BH.
              IF (BODY(I1).GT.0.5*BODY(I)) THEN
                  IF (I1.EQ.2*IPAIR-1) THEN
                      I1 = 2*IPAIR
                  ELSE
                      I1 = 2*IPAIR - 1
                  END IF
              END IF
*       Terminate KS binary and assign kick velocity to single star #I.
              I = I1 + 2*(NPAIRS - IPAIR)
      IF (BODY(I1).GT.0.5*BODY(N+IPAIR)) STOP
              CALL KSTERM
              CALL KICK(I,1,KW1,DM)
              write(*,*) "asigning kick to single object"
*     Initialize new KS polynomials after velocity kick (H > 0 is OK).
              ICOMP = IFIRST
              JCOMP = IFIRST + 1
              CALL KSREG
              IPHASE = -1
*
*       Avoid negative radial velocity for same star (TD2 is small).
              IF (ISKIP.GT.0.AND.TDOT2(NPAIRS).LT.0.0) THEN
                  IF (ECC.GT.1.0) TDOT2(NPAIRS) = 0.0
              END IF
*       Provide diagnostic output if binary survives the kick.
              I = NTOT
              KSTAR(I) = 0
              IF (H(NPAIRS).LT.0.0) THEN
                  SEMI = -0.5d0*BODY(I)/H(NPAIRS)
                  WRITE (6,35)  KW1, R(NPAIRS)/SEMI, SEMI*SU,
     &                          BODY(I)*ZMBAR, (XDOT(K,I)*VSTAR,K=1,3)
   35             FORMAT (' WD/NS BINARY    KW R/A A M V ',
     &                                      I4,3F7.2,3F7.1)
              END IF
          ELSE
*       Set new binary indicator and Roche look-up time (ECC > 0 is OK).
              KSTAR(I) = 0
              CALL TRFLOW(IPAIR,DTR)
              TEV(I) = TIME + DTR
              TMDOT = MIN(TEV(I),TMDOT)
*
*       Update KSTAR & chaos variables for SYNCH after circularization.
              IF (ECC.LE.0.001D0) THEN
                 KSTAR(I) = 10
                 NCE = NCE + 1
                 IC = 0
*       Note: NAMEC may be set in routine CHAOS without call to SPIRAL.
                 DO 36 L = 1,NCHAOS
                     IF (NAME(I).EQ.NAMEC(L)) IC = L
   36            CONTINUE
*       Include safety test on no membership just in case (not sure!).
                 IF (IC.EQ.0) THEN
                     NCHAOS = NCHAOS + 1
                     IC = NCHAOS
                     NAMEC(IC) = NAME(I)
                     IF (NCHAOS.GT.NTMAX) THEN
                          WRITE (6,38)  NAME(I1), NCHAOS, ECC
   38                     FORMAT (' FATAL ERROR!    EXPEL    NM NCH E ',
     &                                              I6,I4,F8.4)
                          STOP
                      END IF
                 END IF
                 IF (IC.GT.0) THEN
                     RP(IC) = SEMI*(1.0 - ECC)
                     ES(IC) = ECC
                     TOSC(IC) = TIME
                 END IF
              END IF
*
*       Include prediction to end of current block-time before new KS.
              TIME = TBLOCK
              IF (T0(I).LT.TIME.AND.DM.GT.0.0D0) THEN
                  CALL XVPRED(I,-2)
                  T0(I) = TIME
                  CALL DTCHCK(TIME,STEP(I),DTK(40))
                  DO 40 K = 1,3
                      X0(K,I) = X(K,I)
                      X0DOT(K,I) = XDOT(K,I)
   40             CONTINUE
              END IF
*
*       Set IPHASE < 0 and re-initialize KS solution with new perturber list.
              IPHASE = -1
              CALL KSLIST(IPAIR)
              CALL KSPOLY(IPAIR,1)
*       Avoid negative radial velocity at pericentre (bug fix 12/08).
              IF (TDOT2(IPAIR).LT.0.0D0) TDOT2(IPAIR) = 0.0
*
*       Correct neighbour forces and update system energy loss due to DM.
              IF (DM.GT.0.0) THEN
                  IF (DM*SMU.LT.0.05) THEN
                      CALL FICORR(I,DM)
                  ELSE
                      CALL FCORR(I,DM,0)
                  END IF
              END IF
*
              IF (DM*SMU.GT.0.1) THEN
*       Include body #I at the end (counting from location #2).
                  NNB2 = NNB + 2
                  ILIST(NNB2) = I
*
*       Obtain new F & FDOT and time-steps.
                  DO 50 L = 2,NNB2
                      J = ILIST(L)
                      IF (L.EQ.NNB2) THEN
                          J = I
                      ELSE
                          CALL XVPRED(J,-2)
                          CALL DTCHCK(TIME,STEP(J),DTK(40))
                          DO 45 K = 1,3
                              X0DOT(K,J) = XDOT(K,J)
   45                     CONTINUE
                      END IF
                      CALL FPOLY1(J,J,0)
                      CALL FPOLY2(J,J,0)
   50             CONTINUE
              END IF
              TPREV = TIME - STEPX
          END IF
      END IF
*
*       Include cleaning of CHAOS arrays (cf. KSTERM for current COAL).
      L = 0
      IC = 0
   52 L = L + 1
      KK = 0
      DO 55 I = N+1,NTOT
          IF (NAME(I).EQ.NAMEC(L)) KK = 1
   55 CONTINUE
      NCH1 = NCHAOS
*       Perform removal of NAMEC on failed search.
      IF (KK.EQ.0.AND.KSTAR(N+IPAIR).LT.0) THEN
          II = -L
          CALL SPIRAL(II)
      END IF
      IC = IC + 1
*       Consider the same location again after array update.
      IF (NCHAOS.LT.NCH1) L = L - 1
      IF (IC.LT.NCHAOS+1) GO TO 52
*
*       Reverse case indicator for exit from collision routine.
   60 ICASE = -ICASE
*
      RETURN
*
      END
