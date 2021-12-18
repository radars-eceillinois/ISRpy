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
        self.n0 = a_WGS / np.sqrt(1 - flatness * (2 - flatness) * np.sin(self.lat0)**2.)
        # cartesian geocentric coordinates wrt Greenwich
        self.x0 = (self.n0 + self.h0) * np.cos(self.lat0) * np.cos(self.lon0)
        self.y0 = (self.n0 + self.h0) * np.cos(self.lat0) * np.sin(self.lon0)
        self.z0 = (self.n0 * (1 - eccentricity**2) + self.h0) * np.sin(self.lat0)
        self.xyz0 = np.array([self.x0, self.y0, self.z0])
        xy0 = np.array([self.x0, self.y0])
        r0 = np.sqrt(np.dot(self.xyz0, self.xyz0))
        p0 = np.sqrt(np.dot(xy0, xy0))

        # unit vectors
        self.east0 = np.array([-self.y0, self.x0, 0]) / p0 # zenith and north directions wrt local ellipsoid
        self.zenith0 = np.array([np.cos(self.lat0) * np.cos(self.lon0),
                                 np.cos(self.lat0) * np.sin(self.lon0),
                                 np.sin(self.lat0)])
        self.north0 = np.cross(self.zenith0, self.east0)

        self.uo = np.array([np.cos(self.dec) * np.cos(self.ha/4. + self.lon0),
                            np.cos(self.dec) * np.sin(self.ha/4. + self.lon0),
                            np.sin(self.dec)])    # on axis
