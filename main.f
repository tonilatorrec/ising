      program main

      use random !generates random integers and real numbers
      use thermo !computes thermodynamic variables

C     ******************************************************************      
      IMPLICIT NONE
C     ******************************************************************
C     VARIABLES

      INTEGER, parameter :: L=30 !lattice size
      INTEGER :: N = L*L !number of spins

      INTEGER SEED, SEED0, NSEED
      DOUBLE PRECISION TEMP, TEMP0, TEMP1, TEMPSTEP
      INTEGER itemp, NTEMP
      INTEGER I, J !control variables

      DOUBLE PRECISION X, MAG, MAGNE, E, ENERG

c     files
      character(len=9) :: outfile ! output file
      character(len=10) :: configfile = "params.nml" 
c-----------------------------------------------------------------------
C     metropolis variables      
      INTEGER I0, J0 !chosen spin
      DOUBLE PRECISION DELTA
      INTEGER DELTAE
      INTEGER MCTOT, MCINI !number of metropolis steps
      INTEGER NITER, IMC, MCD
      DOUBLE PRECISION W(-8:8)

      integer(kind=2) :: s(l,l) !spin matrix
      INTEGER PBC(0:L+1)
c-----------------------------------------------------------------------
C     sums and averages
      DOUBLE PRECISION :: sum=0.d0, SUME=0.d0, SUME2=0.d0, SUMM=0.d0,
     &SUMAM=0.d0, SUMM2=0.d0
      double precision :: vare, varm

      double precision :: EPREV=0.d0 !previous energy (to compute dE/dT)
      DOUBLE PRECISION SUMEN,SUMMN,SUMAMN,VARMN,EPSE,EPSM
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------
C     thermodynamical variables
      DOUBLE PRECISION CAP, SUSC
c-----------------------------------------------------------------------     
      ! Namelist definition.
      namelist /simulparams/ temp0, temp1, ntemp, seed0, nseed
      namelist /mcparams/ mctot, mcini, mcd

      open(20, file=configfile)
      read(20, nml=simulparams)
      read(20, nml=mcparams)
      close(20)

      TEMPSTEP = (TEMP1-TEMP0)/DBLE(NTEMP)
      WRITE(outfile, '(A3,I2,A4)') '_L_',L,'.txt'
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      OPEN(UNIT=13, FILE=outfile, FORM='FORMATTED')
      WRITE(13,104) "TEMP","SUMEN","EPSE","SQRT(SUMM2/N)",
     &"SUMAMN","EPSM","CAP","SUSC", "DE/DT"

c     begin temperature loop
c      DO TEMP=TEMP0,TEMP1,TEMPSTEP
      do itemp=0,ntemp
        temp=temp0+itemp*tempstep
        WRITE(*,'(a, F5.3)') "TEMP = ", TEMP

c     begin seed loop
      DO SEED=SEED0,SEED0+NSEED-1,1
        CALL INIT_GENRAND(SEED)
c-----------------------------------------------------------------------
C     generate initial spin matrix
      DO J=1,L  
      DO I=1,L
          X = GENRAND_REAL2()
          IF (X.LT.0.5D0) THEN
          S(I,J)=1
          ELSE
          S(I,J)=-1
          ENDIF
        ENDDO
        ENDDO
C     specify periodic boundary conditions
        PBC(0)=L
        PBC(L+1)=1
        DO I=1,L
          PBC(I)=I
        ENDDO
C     initial magnetization and energy
        MAG = MAGNE(S,L)
        E = ENERG(S,L,PBC)
c-----------------------------------------------------------------------
C     METROPOLIS ALGORITHM
c-----------------------------------------------------------------------
c     prepare all relevant exponential factors
      DO DELTAE=-8,8
          W(DELTAE)=EXP(-DELTAE/TEMP)
        ENDDO

c     begin metropolis loop
      DO IMC=1,MCTOT
          DO NITER=1,N
            I = 1
            J = L
c           pick random spin at position (i0, j0)
            CALL RANDINT(I,J,I0)
            CALL RANDINT(I,J,J0)
c           compute change in energy if the spin was flipped            
            DELTAE = 2*S(I0,J0)*(S(I0,PBC(J0+1))+S(I0,PBC(J0-1))+
     & S(PBC(I0+1),J0)+S(PBC(I0-1),J0))
            IF (DELTAE.LE.0.D0) THEN ! flip if energy decreases
              S(I0,J0)=-S(I0,J0) 
              E = E+DELTAE
              IF (S(I0,J0).EQ.1) THEN
                MAG = MAG+2
              ELSE
                MAG = MAG-2
              ENDIF
            ELSE ! flip with probability exp(-deltaE/temp)
              DELTA = GENRAND_REAL2()
              IF (DELTA.LT.W(DELTAE)) THEN
                S(I0,J0)=-S(I0,J0)
                E = E+DELTAE
                IF (S(I0,J0).EQ.1) THEN
                  MAG = MAG+2
                ELSE
                  MAG = MAG-2
                ENDIF
              ENDIF
            ENDIF
          ENDDO

c     compute averages every number of steps
          IF ((IMC.GT.MCINI).AND.(MCD*(IMC/MCD).EQ.IMC)) THEN
            MAG = MAGNE(S,L)
            SUM=SUM+1
            SUME=SUME+E
            SUME2=SUME2+E*E
            SUMM=SUMM+MAG
            SUMAM=SUMAM+ABS(MAG)
            SUMM2=SUMM2+MAG*MAG
          ENDIF      
        ENDDO
c     end metropolis loop
c-----------------------------------------------------------------------
      ENDDO
c     end seed loop
c-----------------------------------------------------------------------      
      SUME=SUME/SUM
      SUME2=SUME2/SUM
      SUMM = SUMM/SUM
      SUMAM=SUMAM/SUM
      SUMM2=SUMM2/SUM
      VARE = SUME2-SUME*SUME
      VARM = SUMM2-SUMAM*SUMAM

      SUMEN = SUME/N
      SUMMN = SUMM/N
      SUMAMN = SUMAM/N
      VARMN = VARM/N
      EPSE = VARE/(N*SQRT(SUM))
      EPSM = VARMN/(N*SQRT(SUM))
      CAP = (SUME2-SUME**2)/(N*TEMP**2)
      SUSC = (SUMM2-SUMAM**2)/(N*TEMP)

C     OUTPUT
      WRITE(13,103) TEMP,SUMEN,EPSE,DSQRT(SUMM2)/N,SUMAMN,EPSM,CAP,
     &SUSC,(SUMEN-EPREV)/(TEMPSTEP)
      WRITE(*,*) 
      EPREV=SUMEN
      ENDDO
c     end temperature loop

c-----------------------------------------------------------------------
C     ******************************************************************
C     FORMATS
C     ******************************************************************
  103 FORMAT(9(F13.4,2X))
  104 FORMAT(9(A13,2X))
C     ******************************************************************
      END PROGRAM

C     ******************************************************************