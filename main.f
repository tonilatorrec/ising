      program main

C     ******************************************************************      
      IMPLICIT NONE
C     ******************************************************************
C     VARIABLES
      INTEGER SEED, SEED0, NSEED
      DOUBLE PRECISION TEMP, TEMP0, TEMP1, TEMPSTEP
      INTEGER NTEMP
      INTEGER I, J, K !control variables
      INTEGER L !lattice size
      PARAMETER(L=30)
      INTEGER N !number of spins
      INTEGER*2 S(1:L,1:L) !spin matrix
      INTEGER*4 PBC(0:L+1)
      DOUBLE PRECISION X, GENRAND_REAL2, MAG, MAGNE, E, ENERG
      CHARACTER*13 FILENAME
c-----------------------------------------------------------------------
C     metropolis variables      
      INTEGER I0, J0 !chosen spin
      DOUBLE PRECISION DELTA
      INTEGER*2 DELTAE
      INTEGER MCTOT, MCINI !number of metropolis steps
      INTEGER NITER, IMC, MCD
      DOUBLE PRECISION W(-8:8)
c-----------------------------------------------------------------------
C     sums and averages
      DOUBLE PRECISION SUM, SUME, SUME2, SUMM, SUMAM, SUMM2, VARE, VARM,
     & EPREV ! previous energy (used to compute dE/dT)
      DOUBLE PRECISION SUMEN,SUMMN,SUMAMN,VARMN,EPSE,EPSM
c-----------------------------------------------------------------------     
c-----------------------------------------------------------------------
C     thermodynamical variables
      DOUBLE PRECISION CAP, SUSC
c-----------------------------------------------------------------------     
      N = L*L
      TEMP0 = 1.4D0
      TEMP1 = 3.4D0
      NTEMP = 10
      TEMPSTEP = (TEMP1-TEMP0)/DBLE(NTEMP)

      SEED0=117654
      NSEED=10

      MCTOT=4000
      MCINI=200
      MCD=100
      EPREV=0.D0

      WRITE(FILENAME, '(A7,I2,A4)') 'res/_L_',L,'.txt'
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      OPEN(UNIT=13, FILE=FILENAME, FORM='FORMATTED')
      WRITE(13,104) "TEMP","SUMEN","EPSE","SQRT(SUMM2/N)",
     &"SUMAMN","EPSM","CAP","SUSC", "DE/DT"

c     begin temperature loop
      DO TEMP=TEMP0,TEMP1,TEMPSTEP
        WRITE(*,*) "TEMP = ", TEMP

C     initialize sums of variables to study
      SUM=0.0D0
      SUME=0.0D0
      SUME2=0.0D0
      SUMM=0.0D0
      SUMM2=0.0D0
      SUMAM=0.0D0

c     begin seed loop
      DO SEED=SEED0,SEED0+NSEED-1,1
        CALL INIT_GENRAND(SEED)
c-----------------------------------------------------------------------
C     generate initial spin matrix
        DO I=1,L
        DO J=1,L
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
            CALL RANDINT(I,J,I0)
            CALL RANDINT(I,J,J0)
            DELTAE = 2*S(I0,J0)*(S(I0,PBC(J0+1))+S(I0,PBC(J0-1))+
     & S(PBC(I0+1),J0)+S(PBC(I0-1),J0))
            IF (DELTAE.LE.0.D0) THEN
              S(I0,J0)=-S(I0,J0)
              E = E+DELTAE
              IF (S(I0,J0).EQ.1) THEN
                MAG = MAG+2
              ELSE
                MAG = MAG-2
              ENDIF
            ELSE
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
c-----------------------------------------------------------------------
c     initialize mt(0:N-1) with a seed
c-----------------------------------------------------------------------
      subroutine init_genrand(s)
      integer s
      integer N
      integer DONE
      integer ALLBIT_MASK
      parameter (N=624)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask1/ ALLBIT_MASK
c
      call mt_initln
      mt(0)=iand(s,ALLBIT_MASK)
      do 100 mti=1,N-1
        mt(mti)=1812433253*
     &          ieor(mt(mti-1),ishft(mt(mti-1),-30))+mti
        mt(mti)=iand(mt(mti),ALLBIT_MASK)
  100 continue
      initialized=DONE
c
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,1)-real-interval
c-----------------------------------------------------------------------
      function genrand_real2()
      double precision genrand_real2,r
      integer genrand_int32
      r=dble(genrand_int32())
      if(r.lt.0.d0)r=r+2.d0**32
      genrand_real2=r/4294967296.d0
      return
      end
c-----------------------------------------------------------------------
c     generates a random number on [0,0xffffffff]-interval
c-----------------------------------------------------------------------
      function genrand_int32()
      integer genrand_int32
      integer N,M
      integer DONE
      integer UPPER_MASK,LOWER_MASK,MATRIX_A
      integer T1_MASK,T2_MASK
      parameter (N=624)
      parameter (M=397)
      parameter (DONE=123456789)
      integer mti,initialized
      integer mt(0:N-1)
      integer y,kk
      integer mag01(0:1)
      common /mt_state1/ mti,initialized
      common /mt_state2/ mt
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
c
      if(initialized.ne.DONE)then
        call init_genrand(21641)
      endif
c
      if(mti.ge.N)then
        do 100 kk=0,N-M-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+M),ishft(y,-1)),mag01(iand(y,1)))
  100   continue
        do 200 kk=N-M,N-1-1
          y=ior(iand(mt(kk),UPPER_MASK),iand(mt(kk+1),LOWER_MASK))
          mt(kk)=ieor(ieor(mt(kk+(M-N)),ishft(y,-1)),mag01(iand(y,1)))
  200   continue
        y=ior(iand(mt(N-1),UPPER_MASK),iand(mt(0),LOWER_MASK))
        mt(kk)=ieor(ieor(mt(M-1),ishft(y,-1)),mag01(iand(y,1)))
        mti=0
      endif
c
      y=mt(mti)
      mti=mti+1
c
      y=ieor(y,ishft(y,-11))
      y=ieor(y,iand(ishft(y,7),T1_MASK))
      y=ieor(y,iand(ishft(y,15),T2_MASK))
      y=ieor(y,ishft(y,-18))
c
      genrand_int32=y
      return
      end
c-----------------------------------------------------------------------
c     initialize large number (over 32-bit constant number)
c-----------------------------------------------------------------------
      subroutine mt_initln
      integer ALLBIT_MASK
      integer TOPBIT_MASK
      integer UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      integer mag01(0:1)
      common /mt_mask1/ ALLBIT_MASK
      common /mt_mask2/ TOPBIT_MASK
      common /mt_mask3/ UPPER_MASK,LOWER_MASK,MATRIX_A,T1_MASK,T2_MASK
      common /mt_mag01/ mag01
CC    TOPBIT_MASK = Z'80000000'
CC    ALLBIT_MASK = Z'ffffffff'
CC    UPPER_MASK  = Z'80000000'
CC    LOWER_MASK  = Z'7fffffff'
CC    MATRIX_A    = Z'9908b0df'
CC    T1_MASK     = Z'9d2c5680'
CC    T2_MASK     = Z'efc60000'
      TOPBIT_MASK=1073741824
      TOPBIT_MASK=ishft(TOPBIT_MASK,1)
      ALLBIT_MASK=2147483647
      ALLBIT_MASK=ior(ALLBIT_MASK,TOPBIT_MASK)
      UPPER_MASK=TOPBIT_MASK
      LOWER_MASK=2147483647
      MATRIX_A=419999967
      MATRIX_A=ior(MATRIX_A,TOPBIT_MASK)
      T1_MASK=489444992
      T1_MASK=ior(T1_MASK,TOPBIT_MASK)
      T2_MASK=1875247104
      T2_MASK=ior(T2_MASK,TOPBIT_MASK)
      mag01(0)=0
      mag01(1)=MATRIX_A
      return
      end
C     ******************************************************************
      DOUBLE PRECISION FUNCTION MAGNE(S,L)
C     ******************************************************************
C     computes system magnetization
C     input:
C     S (INTEGER*2): spin matrix 
C     L (INTEGER*4): matrix size
C     ******************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION MAG
      INTEGER L,I,J
      INTEGER*2 S(1:L,1:L)
C     ******************************************************************
      MAG = 0.D0
      DO J=1,L
      DO I=1,L
        MAG = MAG+S(I,J)
      ENDDO
      ENDDO
      MAGNE = MAG
      RETURN
      END
C     ******************************************************************
      DOUBLE PRECISION FUNCTION ENERG(S,L,PBC)
C     ******************************************************************
C     computes system energy
C     input:
C     S (INTEGER*2): spin matrix 
C     L (INTEGER*4): matrix size
C     PBC (INTEGER*4): boundary conditions matrix
C     ******************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION E
      INTEGER L,I,J,INT
      INTEGER*2 S(1:L,1:L)
      INTEGER*4 PBC(0:L+1)
C     ******************************************************************
      E=0.D0

      DO J=1,L
      DO I=1,L
        INT = -S(PBC(I-1),J)*S(I,J)-S(I,PBC(J-1))*S(I,J)
        E = E+INT
      ENDDO
      ENDDO

      ENERG = E
      RETURN
      END
c-----------------------------------------------------------------------
      SUBROUTINE RANDINT(A,B,N)
c----------------------------------------------------------------------- 
      IMPLICIT NONE
      INTEGER*4 A,B,N
      DOUBLE PRECISION X, GENRAND_REAL2
c-----------------------------------------------------------------------
      X = GENRAND_REAL2()
      N = FLOOR(X*(1+B-A))+A
      RETURN
      END

