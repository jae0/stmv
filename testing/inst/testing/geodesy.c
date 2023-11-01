/* geodesy.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Subroutine */ int inver1_(doublereal *dlat1, doublereal *dlon1, doublereal 
	*dlat2, doublereal *dlon2, doublereal *a, doublereal *f, doublereal *
	faz, doublereal *baz, doublereal *s)
{
    /* Initialized data */

    static doublereal eps = 5e-14;
    static doublereal pi = 3.1415926535897932384626433832795;
    static doublereal rad = .0174532925199432957692369076848861;

    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sin(doublereal), cos(doublereal), sqrt(doublereal), atan2(
	    doublereal, doublereal);

    /* Local variables */
    static doublereal c__, d__, e, r__, x, y, sa, cx, cy, cz, sx, sy, c2a, 
	    cu1, cu2, su1, tu1, tu2, glat1, glat2, glon1, glon2;


/* *** SOLUTION OF THE GEODETIC INVERSE PROBLEM AFTER T.VINCENTY */
/* *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS */
/* *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL */
/* *** STANDPOINT/FOREPOINT MUST NOT BE THE GEOGRAPHIC POLE */

/* *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID */
/* *** F IS THE FLATTENING (NOT RECIPROCAL) OF THE REFERNECE ELLIPSOID */
/* *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST */
/* *** FORWARD AZIMUTHS AT BOTH POINTS RETURNED IN RADIANS FROM NORTH */

/* *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 18FEB75 */
/* *** MODIFIED FOR IBM SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 7507 */

/* *** Modified for R by D.Gillis Zoology University of Manitoba 16JUN03 */
/*     Replaced common blocks for constants with DATA statements.  Datum */
/*     parameters moved from common block to subroutine arguements. */

/* *** Input Variables: DLAT1,DLON1 - initial fix (P1) in degrees */
/*                          (latitude, north +) (longitude east +) */
/*                      DLAT2,DLON2 - destination fix (P2) in degrees */

/*                      Ellipsoid (spheroid model, eg. WGS84) */
/*                         A   - radius of major axis in distance units */
/*                         F   - flattening factor */

/* *** Output Variables: FAZ - azimuth of the geodesic (P1 to P2) */
/*                       BAZ - azimuth of the geodesic (P2 to P1) */
/*                       S   - spheroidal distance = length of the geodesic */
/*                             in distance units */

/* *** After Vincenty,T. 1975. Direct and inverse solutions of geodesics */
/*           on the ellipsoid with application of nested equations. Survey */
/*           Review 23(176):88-94. */


/*     These routines were compiled in Windows 98 (DOS Window) using the */
/*     Gnu FORTRAN complier: g77 --share -o geodesy.dll marspat.for */
/*     creating:             geodesy.dll */



/*     COMMON/CONST/PI,RAD - original code */
/*     COMMON/ELIPSOID/A,F - original code */
    r__ = 1. - *f;

/*     Convert decimal degrees to radians */

    glat1 = *dlat1 * rad;
    glon1 = *dlon1 * rad;
    glat2 = *dlat2 * rad;
    glon2 = *dlon2 * rad;

    tu1 = r__ * sin(glat1) / cos(glat1);
    tu2 = r__ * sin(glat2) / cos(glat2);
    cu1 = 1.f / sqrt(tu1 * tu1 + 1.f);
    su1 = cu1 * tu1;
    cu2 = 1.f / sqrt(tu2 * tu2 + 1.f);
    *s = cu1 * cu2;
    *baz = *s * tu2;
    *faz = *baz * tu1;
    x = glon2 - glon1;
L100:
    sx = sin(x);
    cx = cos(x);
    tu1 = cu2 * sx;
    tu2 = *baz - su1 * cu2 * cx;
    sy = sqrt(tu1 * tu1 + tu2 * tu2);
    cy = *s * cx + *faz;
    y = atan2(sy, cy);
    sa = *s * sx / sy;
    c2a = -sa * sa + 1.f;
    cz = *faz + *faz;
    if (c2a > 0.f) {
	cz = -cz / c2a + cy;
    }
    e = cz * cz * 2.f - 1.f;
    c__ = ((c2a * -3.f + 4.f) * *f + 4.f) * c2a * *f / 16.f;
    d__ = x;
    x = ((e * cy * c__ + cz) * sy * c__ + y) * sa;
    x = (1.f - c__) * x * *f + glon2 - glon1;
    if ((d__1 = d__ - x, abs(d__1)) > eps) {
	goto L100;
    }
    *faz = atan2(tu1, tu2);
    *baz = atan2(cu1 * sx, *baz * cx - su1 * cu2) + pi;
    x = sqrt((1.f / r__ / r__ - 1.f) * c2a + 1.f) + 1.f;
    x = (x - 2.f) / x;
    c__ = 1.f - x;
    c__ = (x * x / 4.f + 1.f) / c__;
    d__ = (x * .375f * x - 1.f) * x;
    x = e * cy;
    *s = 1.f - e - e;
    *s = ((((sy * sy * 4.f - 3.f) * *s * cz * d__ / 6.f - x) * d__ / 4.f + cz)
	     * sy * d__ + y) * c__ * *a * r__;
    *faz /= rad;
    *baz /= rad;
    return 0;
} /* inver1_ */

