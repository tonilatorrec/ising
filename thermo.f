      module thermo

      implicit none

      contains

      subroutine thermo_params(s, l, mode)
c     ------------------------------------------------------------------

c     ------------------------------------------------------------------
      implicit none
      integer, parameter :: short = selected_int_kind(2)
     
      integer :: l
      integer :: mode !1 for magnetization, 2 for energy
      integer :: i, j
      integer(kind=2) s

      double precision :: mag = 0.d0, e = 0.d0
c     ------------------------------------------------------------------

      end

      DOUBLE PRECISION FUNCTION MAGNE(S,L)
C     ******************************************************************
C     computes system magnetization
C     input:
C     S (INTEGER): spin matrix 
C     L (INTEGER): matrix size
C     ******************************************************************
      IMPLICIT NONE
      integer, parameter :: short = selected_int_kind(2)

      DOUBLE PRECISION MAG
      INTEGER L,I,J
      integer(kind=2) :: S(L,L)
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
C     S (INTEGER, short): spin matrix 
C     L (INTEGER): matrix size
C     PBC (INTEGER): boundary conditions matrix
C     ******************************************************************
      IMPLICIT NONE
      integer, parameter :: short = selected_int_kind(2)

      DOUBLE PRECISION E
      INTEGER L,I,J,INT
      integer(kind=2) s(l,l)
      INTEGER PBC(0:L+1)
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

      end module