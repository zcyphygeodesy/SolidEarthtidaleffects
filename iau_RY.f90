      SUBROUTINE iau_RY ( THETA, R )
!+
!  - - - - - - -
!   i a u _ R Y
!  - - - - - - -
!
!  Rotate an r-matrix about the y-axis.
!
!  This routine is part of the International Astronomical Union's
!  SOFA (Standards of Fundamental Astronomy) software collection.
!
!  Status:  vector/matrix support routine.
!
!  Given:
!     THETA    d         angle (radians)
!
!  Given and returned:
!     R        d(3,3)    r-matrix
!
!  Sign convention:  The matrix can be used to rotate the
!  reference frame of a vector.  Calling this routine with
!  positive THETA incorporates in the matrix an additional
!  rotation, about the y-axis, anticlockwise as seen looking
!  towards the origin from positive y.
!
!  Called:
!     iau_IR       initialize r-matrix to identity
!     iau_RXR      r-matrix multiply
!     iau_CR       r-matrix copy
!
!  This revision:  2000 November 25
!
!  Copyright (C) 2003 IAU SOFA Review Board.  See notes at end.
!
!-----------------------------------------------------------------------

      IMPLICIT NONE

      REAL*8 THETA, R(3,3)

      REAL*8 S, C, A(3,3), W(3,3)

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

!  Matrix representing new rotation.
      S = SIN(THETA)
      C = COS(THETA)
      CALL iau_IR ( A )
      A(1,1) = C
      A(3,1) = S
      A(1,3) = -S
      A(3,3) = C

!  Rotate.
      CALL iau_RXR ( A, R, W )

!  Return result.
      CALL iau_CR ( W, R )

!  Finished.

!+----------------------------------------------------------------------
!
!  Copyright (C) 2003
!  Standards Of Fundamental Astronomy Review Board
!  of the International Astronomical Union.
!
!  =====================
!  SOFA Software License
!  =====================
!
!  NOTICE TO USER:
!
!  BY USING THIS SOFTWARE YOU ACCEPT THE FOLLOWING TERMS AND CONDITIONS
!  WHICH APPLY TO ITS USE.
!
!  1. The Software is owned by the IAU SOFA Review Board ("the Board").
!
!  2. The Software is made available free of charge for use by:
!
!     a) private individuals for non-profit research; and
!
!     b) non-profit educational, academic and research institutions.
!
!  3. Commercial use of the Software is specifically excluded from the
!     terms and conditions of this license.  Commercial use of the
!     Software is subject to the prior written agreement of the Board on
!     terms to be agreed.
!
!  4. The provision of any version of the Software under the terms and
!     conditions specified herein does not imply that future versions
!     will also be made available under the same terms and conditions.
!
!  5. The user may modify the Software for his/her own purposes.  The
!     user may distribute the modified software provided that the Board
!     is informed and that a copy of the modified software is made
!     available to the Board on request.  All modifications made by the
!     user shall be clearly identified to show how the modified software
!     differs from the original Software, and the name(s) of the
!     affected routine(s) shall be changed.  The original SOFA Software
!     License text must be present.
!
!  6. In any published work produced by the user and which includes
!     results achieved by using the Software, the user shall acknowledge
!     that the Software was used in producing the information contained
!     in such publication.
!
!  7. The user may incorporate or embed the Software into other software
!     products which he/she may then give away free of charge but not
!     sell provided the user makes due acknowledgement of the use which
!     he/she has made of the Software in creating such software
!     products.  Any redistribution of the Software in this way shall be
!     made under the same terms and conditions under which the user
!     received it from the SOFA Center.
!
!  8. The user shall not cause the Software to be brought into
!     disrepute, either by misuse, or use for inappropriate tasks, or by
!     inappropriate modification.
!
!  9. The Software is provided to the user "as is" and the Board makes
!     no warranty as to its use or performance.   The Board does not and
!     cannot warrant the performance or results which the user may
!     obtain by using the Software.  The Board makes no warranties,
!     express or implied, as to non-infringement of third party rights,
!     merchantability, or fitness for any particular purpose.  In no
!     event will the Board be liable to the user for any consequential,
!     incidental, or special damages, including any lost profits or lost
!     savings, even if a Board representative has been advised of such
!     damages, or for any claim by any third party.
!
!  Correspondence concerning SOFA software should be addressed as
!  follows:
!
!     Internet email: sofa@rl.ac.uk
!     Postal address: IAU SOFA Center
!                     Rutherford Appleton Laboratory
!                     Chilton, Didcot, Oxon OX11 0QX
!                     United Kingdom
!
!
!-----------------------------------------------------------------------

      END
