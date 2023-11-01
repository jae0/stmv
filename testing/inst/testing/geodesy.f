      SUBROUTINE INVER1(DLAT1,DLON1,DLAT2,DLON2,A,F,FAZ,BAZ,S)
C
C *** SOLUTION OF THE GEODETIC INVERSE PROBLEM AFTER T.VINCENTY
C *** MODIFIED RAINSFORD'S METHOD WITH HELMERT'S ELLIPTICAL TERMS
C *** EFFECTIVE IN ANY AZIMUTH AND AT ANY DISTANCE SHORT OF ANTIPODAL
C *** STANDPOINT/FOREPOINT MUST NOT BE THE GEOGRAPHIC POLE
C                              
C *** A IS THE SEMI-MAJOR AXIS OF THE REFERENCE ELLIPSOID
C *** F IS THE FLATTENING (NOT RECIPROCAL) OF THE REFERNECE ELLIPSOID
C *** LATITUDES AND LONGITUDES IN RADIANS POSITIVE NORTH AND EAST
C *** FORWARD AZIMUTHS AT BOTH POINTS RETURNED IN RADIANS FROM NORTH
C
C *** PROGRAMMED FOR CDC-6600 BY LCDR L.PFEIFER NGS ROCKVILLE MD 18FEB75
C *** MODIFIED FOR IBM SYSTEM 360 BY JOHN G GERGEN NGS ROCKVILLE MD 7507
C
C *** Modified for R by D.Gillis Zoology University of Manitoba 16JUN03
C     Replaced common blocks for constants with DATA statements.  Datum
C     parameters moved from common block to subroutine arguements.
C
C *** Input Variables: DLAT1,DLON1 - initial fix (P1) in degrees
C                          (latitude, north +) (longitude east +)
C                      DLAT2,DLON2 - destination fix (P2) in degrees
C
C                      Ellipsoid (spheroid model, eg. WGS84)
C                         A   - radius of major axis in distance units
C                         F   - flattening factor
C
C *** Output Variables: FAZ - azimuth of the geodesic (P1 to P2)
C                       BAZ - azimuth of the geodesic (P2 to P1)
C                       S   - spheroidal distance = length of the geodesic
C                             in distance units
C
C *** After Vincenty,T. 1975. Direct and inverse solutions of geodesics
C           on the ellipsoid with application of nested equations. Survey
C           Review 23(176):88-94.
C
C
C     These routines were compiled in Windows 98 (DOS Window) using the 
C     Gnu FORTRAN complier: g77 --share -o geodesy.dll marspat.for
C     creating:             geodesy.dll
C
C
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C     COMMON/CONST/PI,RAD - original code
C     COMMON/ELIPSOID/A,F - original code
      DATA EPS/0.5D-13/
      DATA PI/3.1415926535897932384626433832795D0/
      DATA RAD/0.0174532925199432957692369076848861D0/
      R=1.D0-F
C
C     Convert decimal degrees to radians
C
      GLAT1=DLAT1*RAD
      GLON1=DLON1*RAD
      GLAT2=DLAT2*RAD
      GLON2=DLON2*RAD
C
      TU1=R*DSIN(GLAT1)/DCOS(GLAT1)
      TU2=R*DSIN(GLAT2)/DCOS(GLAT2)
      CU1=1./DSQRT(TU1*TU1+1.)
      SU1=CU1*TU1
      CU2=1./DSQRT(TU2*TU2+1.)
      S=CU1*CU2
      BAZ=S*TU2
      FAZ=BAZ*TU1
      X=GLON2-GLON1
  100 SX=DSIN(X)
      CX=DCOS(X)
      TU1=CU2*SX
      TU2=BAZ-SU1*CU2*CX
      SY=DSQRT(TU1*TU1+TU2*TU2)
      CY=S*CX+FAZ
      Y=DATAN2(SY,CY)
      SA=S*SX/SY
      C2A=-SA*SA+1.
      CZ=FAZ+FAZ
      IF(C2A.GT.0.)CZ=-CZ/C2A+CY
      E=CZ*CZ*2.-1.
      C=((-3.*C2A+4.)*F+4.)*C2A*F/16.
      D=X
      X=((E*CY*C+CZ)*SY*C+Y)*SA
      X=(1.-C)*X*F+GLON2-GLON1
      IF(DABS(D-X).GT.EPS) GOTO 100
      FAZ=DATAN2(TU1,TU2)
      BAZ=DATAN2(CU1*SX,BAZ*CX-SU1*CU2)+PI
      X=DSQRT((1./R/R-1.)*C2A+1.)+1.
      X=(X-2.)/X
      C=1.-X
      C=(X*X/4.+1.)/C
      D=(0.375*X*X-1.)*X
      X=E*CY
      S=1.-E-E
      S=((((SY*SY*4.-3.)*S*CZ*D/6.-X)*D/4.+CZ)*SY*D+Y)*C*A*R
      FAZ=FAZ/RAD
      BAZ=BAZ/RAD
      END

