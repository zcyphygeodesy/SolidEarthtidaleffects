## Fortran codes for computation of solid Earth tidal effects on all-element geodetic variations on or outside geoid
https://www.zcyphygeodesy.com/en/h-nd-117.html
## [Algorithm purpose]
    Given the longitude, latitude, height and forecast time of the calculation point, compute the solid tidal effects on the geoid or height anomaly (mm), ground gravity (μGal), gravity disturbance (μGal), ground tilt (SW, to the south and to the west, mas), vertical deflection (SW, to the south and to the west, mas), horizontal displacement (EN, to the east and to the north, mm), ground radial displacement (mm), ground normal or orthometric height (mm), radial gravity gradient (10μE) or horizontal gravity gradient (NW, to the north and to the west, 10μE).
    The algorithm adopts the same numerical standards and analytical algorithms of the solid Earth tidal effects on geopotential and geodetic site displacement that are compatible with the IERS Conventions (2010), and compute uniformly the solid Earth tidal effects on all-element geodetic variations considering the latitude and frequency dependence of the Love numbers, so as to maintain rigorously the analytical relationships between the solid Earth tidal effects on various geodetic variations.
    The Earth's tide generating potential (TGP) from the moon is calculated from 2nd to 6th degree, that from the sun from 2nd to 3rd degree, and that from other planets at the 2nd degree.
## [Computation Output]
    tdn(14): the solid tidal effects on all-element geodetic variations.
    tdn(1:14) stores the solid tidal effects on 10 kinds of geodetic variations, which are the solid tidal effects on height anomaly tdn(1) (mm), ground gravity #tdn(2) (μGal), gravity disturbance tdn(3) (μGal), ground tilt #tdn(4:5) (SW, to the south and to the west, mas), vertical deflection tdn(6:7) (SW, to the south and to the west, mas), horizontal displacement #tdn(8:9) (EN, to the east and to the north, mm), ground radial displacement #tdn(10) (mm), ground normal or orthometric height #tdn(11) (mm), radial gravity gradient tdn(12 )(10μE) or horizontal gravity gradient tdn(13:14) (NW, to the north and to the west, 10μE).
    The calculation point can be on the ground, low altitude, satellite, ocean or underwater space. The geodetic variations abvove marked with # are valid only when the site is fixed with the solid Earth.
## [Geophysical models]
    (1) The JPL Planetary Ephemeris DE440 file JEPH.440 (from 1850 to 2201).
    (2) The IERS Earth orientation parameters (EOP) time series file IERSeopc04.dat.
    (3) The Love number correction file frqadjlovekhl.txt for frequency dependence. was generated from Table 6.5a, 6.5b, 6.5c, 7.2, 7.3a and 7.3b in IERS conventions (2010) as well as the maximum amplitude of the equilibrium tidal height of the tidal constituent.
## [Main program for test entrance]
    Solidtidaleffects.f90
    The record of the test output file reslt.txt: the long integer time agreed by ETideLoad, difference between the MJD day and starting MJD0, tdn(1:14)
## (1) Algorithm module for the solid tidal effects on all-element geodetic variations.
    GravFdSldTd(td, BLH, tdn, eoput1, GRS, val, NUM, frd, TRStCRS)
    Input parameters: td - Julian Date (JD), =mjd+2400000.5+32.184/86400.
    Input parameters: BLH(3) - latitude, longitude (decimal degrees) and ellipsoidal height (m) at the calculation point.
    Input parameters: eoput1(6) – EOP [x(") y(") UT1-UT(s) LOD(s) dX(") dY(")] at epoch time td.
    Input parameters: GRS(6)－gm,ae,j2,omega,1/f, default value.
    Input parameters: NUM = 645 for DE440, =156 for DE405.
    Input parameters: TRStCRS(3,3) - conversion matrix of ITRS to GCRS..
    Return parameters: tdn(1:14) - the solid tidal effects on all-element geodetic variations.
## (2) Algorithm module for the geopotential coefficient variations of the Earth's tidal generating potential (TGP) from celestial bodies
    NormTide(td, val, nn, cnm, snm, rln, GRS, TRStCRS)
    Input parameters: cnm, snm - spherical coordinatesthe direct influence of TGP to geopotential coefficient variations.
    Input parameters: rln(3) – the spherical coordinates of calculation point in IERS.
    Return parameters: cnm, snm- the geopotential coefficient variations of the Earth's tidal generating potential (TGP) from celestial bodies, that is the direct influence of TGP to geopotential coefficient variations.
## (3) Algorithm module for the solid tidal effects on all-element geodetic variations using the nominal Love numbers.
    solidtdnormal(rln, 6, cnm, snm, tdn, GRS)
    Input parameters: cnm, snm - the direct influence of TGP to geopotential coefficient variations.
    Return parameters: tdn(1:14) - the solid tidal effects on all-element geodetic variations using the nominal Love numbers.
## (4) Algorithm module for the geopoential coefficient adjustments caused by Love number frequency dependence
    Freqdploven(td,cnm,snm,frd,eop)
    Input parameters: frd(119,8) – the adjustment parameters of Love number frequency dependence.
    Return parameters: cnm, snm - the geopoential coefficient variations caused by Love number frequency dependence.
## (5) Algorithm module for the solid tidal effect adjustments caused by Love number frequency dependence
    Tdlovefreqadj(rln, cnm, snm, ftd, GRS)
    Input parameters: cnm, snm - the geopoential coefficient variations caused by Love number frequency dependence.
    Return parameters: ftd(1:14)- the solid tidal effect adjustments caused by Love number frequency dependence.
## (6) Algorithm module for spherical harmonic synthesis of the anomlous Earth’s gravity field
    RentGField(maxn, rln, cnm, snm, gvm, GRS)
    Return parameters: gvm(9)- disturbance potential(m2/s2), height anomaly (m), gravity anomaly (mGal), gravity disturbance (mGal), vertical deflection vector (ʺ, south, west), disturbing gravity gradient (E, radial), tangential gravity gradient vector (E, north, west).
## (7) Algorithm library for the JPL Planetary Ephemeris
    ReadDE405.f90    !Modified, suitable for DE405 / DE440.
Please download the file JPLEPH.440 into the current directory https://download.s21i.co99.net/24192633/0/0/ABUIABAAGAAg1J_sugYokL_auQc?f=JPLEPH.440&v=1732972524
## (8) Calculation module for normal Earth’s gravity field
    normdjn(GRS, djn); GNormalfd(BLH, NFD, GRS)
## (9) Algorithm module for normalized associative Legendre function and its derivative to θ
    BelPnmdt(pnm, dpt1, dpt2, maxn, t)  ! t = cosθ
## (10) Calculation module for Legendre functions and their derivatives to ψ
    PlmBar_d(p, dp, lmax, rlat) ! ψ = rlat
## (11) Algorithm library for transforming of geodetic coordinates
    BLH_RLAT(GRS, BLH, RLAT); BLH_XYZ(GRS, BLH, XYZ);
    RLAT_BLH(GRS, RLAT, BLH)
## (12) Algorithm library for converting of time system
    CAL2JD (IY0, IM0, ID0, DJM, J); JD2CAL(DJ1, DJ2, IY, IM, ID, FD, J)
## (13) IAU SOFA2000 library
    TRStoCRS(tm, eoput1, TRStCRS); GMST2000(UTA, UTB, TTA, TTB);
    RBPN2000(X, Y, S, RBPN)；……
## (14) Other auxiliary modules
    PickRecord(str0, kln, rec, nn); tmcnt(tm, iyr, imo, idy, ihr, imn, sec)
    mjdtotm(mjd0, ltm); tmtostr(tm, tmstr)
## [For compile and link]
    Fortran90, 132 Columns fixed format. Fortran compiler for any operating system. No external link library required.
## [Algorithmic formula] ETideLoad4.5 User Reference https://www.zcyphygeodesy.com/en/
    8.18.1 Solid tidal effects on various geodetic variations outside solid Earth
    8.2.3 The normalized associated Legendre functions and thier derivative to θ
    8.3.3 Legendre function and its first and second derivatives to ψ
The zip compression package includes the test project in visual studio 2017 - intel fortran integrated environment, DOS executable test file, geophysical models and all input and output data.
