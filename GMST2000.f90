      REAL*8 FUNCTION GMST2000 ( UTA, UTB, TTA, TTB )
*+
*  - - - - - - - - -
*   G M S T 2 0 0 0
*  - - - - - - - - -
*
*  Greenwich Mean Sidereal Time (model consistent with IAU 2000
*  resolutions).
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Given:
*     UTA, UTB     d      UT1 date (JD = UTA+UTB)
*     TTA, TTB     d      TT date (JD = TTA+TTB)
*
*  The result is the Greenwich Mean Sidereal Time (radians), in the
*  range 0 to 2pi.
*
*  Calls SOFA routine iau_ANP and IERS routine ERA2000.
*
*  This revision:  2002 December 2
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL*8 UTA, UTB, TTA, TTB

*  Arcseconds to radians
      REAL*8 DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000), JD
      REAL*8 DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian century
      REAL*8 DJC
      PARAMETER ( DJC = 36525D0 )

      REAL*8 T

      REAL*8 iau_ANP, ERA2000

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  TT Julian centuries since J2000.0.
      T = ( ( TTA-DJ0 ) + TTB ) / DJC

*  Greenwich Mean Sidereal Time, IAU 2000.
      GMST2000 = iau_ANP ( ERA2000 ( UTA, UTB ) +
     :                        (    0.014506  D0 +
     :                        ( 4612.15739966D0 +
     :                        (  + 1.39667721D0 +
     :                        (  - 0.00009344D0 +
     :                        (  + 0.00001882D0 )
     :                                 * T ) * T ) * T ) * T ) * DAS2R )

      END
