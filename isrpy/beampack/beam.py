# -*- coding: utf-8 -*-
"""beam.py

   module for calculating geometry parameters and magnetic aspect
   angle of radar targets

   use aspect_elaz or aspect_txty to calculate aspect angles of targets
   specified by (el,az) or (tx,ty) angles

  Created by Erhan Kudeki on 11/29/08.
  Copyright (c) 2008 ECE, UIUC. All rights reserved.

  History
  Pablo Reyes on 12/15/2021
    - flattening has been changed from 1/298.257 to 1./298.257223563
      using the WGS84 reference in:
      https://gis.icao.int/eganp/webpdf/REF08-Doc9674.pdf
    - First crack at PEP 8 compliance
    - Using NumPy Style Python Docstrings
"""

import numpy as np

# module level variables:

eps = np.finfo(float).eps                  # float resolution
deg = np.pi / 180.                         # to express angles in degree values
a_igrf = 6371.2                            # mean earth radius (km)

a_WGS = 6378.137                           # equatorial radius WGS 84
flatness = 1/298.257223563
b_WGS = a_WGS * (1 - flatness)             # WGS polar radius
eccentricity = np.sqrt(a_WGS**2 - b_WGS**2) / a_WGS

def llh2xyz(latg, lon, h):
    """Converts a point in geodetic to ECEF coordinates.

    Parameters
    ----------
    latg : float
        The geodetic latitude in radians.
    lon : float
        The longitude in radians.
    h : float
        The height above the local ellipsoid

    Returns
    -------
    x : float
        ECEF x coordinate
    y : float
        ECEF y coordinate
    z : float
        ECEF z coordinate
    """

    n = a_WGS / np.sqrt(1 - flatness * (2 - flatness) * np.sin(latg)**2.)
    # cartesian geocentric coordinates wrt Greenwich
    x = (n + h) * np.cos(latg) * np.cos(lon)
    y = (n + h) * np.cos(latg) * np.sin(lon)
    z = (n * (1 - eccentricity**2.) + h) * np.sin(latg)

    return x, y, z

def xyz2llh(x, y, z):
    """Converts a point in  ECEF to geodetic coordinates.

    Returns longitude 'lon', geodetic latitude 'lat', and height 'h'
    of position (x,y,z) defined in geocentric coordinate system

    Parameters
    ----------
    x : float
        ECEF x coordinate
    y : float
        ECEF y coordinate
    z : float
        ECEF z coordinate

    Returns
    -------
    lat : float
        The geodetic latitude in radians.
    lon : float
        The longitude in radians.
    h : float
        The height above the local ellipsoid
    """

    p = np.sqrt(x**2. + y**2.)
    lon = np.arctan2(y, x)
    lat = np.arctan2(z, p)
    latp = lat
    for i in range(10):
        n = a_WGS / np.sqrt(1 - flatness * (2 - flatness) * np.sin(latp)**2.)
        h = p / np.cos(latp) - n
        lat = np.arctan(z / (p * (1 - n * eccentricity**2. / (n + h))))
        if abs(lat - latp) < 3 * eps:
            n = a_WGS / np.sqrt(1 - flatness * (2 - flatness) * np.sin(lat)**2.)
            h = p / np.cos(lat) - n
            break
        latp = lat
    return lat, lon, h

class RadarSpecs:
    """Will contain radar coordinates and coordinate conversions"""

    def __init__(self,lat0, lon0, h0, dec, ha):

        self.lat0 = lat0
        self.lon0 = lon0
        self.h0 = h0
        self.dec = dec
        self.ha = ha

        self.calculate_geometry()

    def calculate_geometry(self):
        self.n0 = a_WGS / np.sqrt(1 - flatness * (2 - flatness) * np.sin(lat0)**2.)
        # cartesian geocentric coordinates wrt Greenwich
        self.x0=(self.n0 + self.h0) * np.cos(self.lat0) * np.cos(self.lon0)
        y0=(n0+h0)*cos(lat0)*sin(lon0)
        z0=(n0*(1-eccentricity**2)+h0)*sin(lat0)
        xyz0=array([x0,y0,z0])
        xy0=array([x0,y0])
        r0=sqrt(dot(xyz0,xyz0))
        p0=sqrt(dot(xy0,xy0))

        # unit vectors from jro
        east0=array([-y0,x0,0])/p0				# zenith and north directions wrt local ellipsoid
        zenith0=array([cos(lat0)*cos(lon0),cos(lat0)*sin(lon0),sin(lat0)])
        north0=cross(zenith0,east0)

        # orthonormal basis vectors including the jro on-axis direction
        dec=-12.88*deg
        ha=-(4+37./60.)*deg						# on-axis direction at JRO
        uo=array([cos(dec)*cos(ha/4.+lon0),cos(dec)*sin(ha/4.+lon0),sin(dec)])	# on axis
        ux=cross(zenith0,uo)
        ux=ux/sqrt(dot(ux,ux))					# along the building to the right
        uy=cross(uo,ux)							# away from the building into the valley










