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
from ..pyigrf import pyigrf

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

class TargetGeometry:
    """Will contain radar coordinates and coordinate conversions"""

    def __init__(self,lat0, lon0, h0, IGRFmodel=None):

        self.lat0 = lat0
        self.lon0 = lon0
        self.h0 = h0

        # None argument instantiates the latest
        self.igrf = pyigrf(IGRFmodel)
        # igrf to be instantiated once and to be remembered  between calls

        n0 = a_WGS / np.sqrt(1 - flatness * (2 - flatness) * np.sin(self.lat0)**2.)
        # cartesian geocentric coordinates wrt Greenwich
        x0 = (n0 + self.h0) * np.cos(self.lat0) * np.cos(self.lon0)
        y0 = (n0 + self.h0) * np.cos(self.lat0) * np.sin(self.lon0)
        z0 = (n0 * (1 - eccentricity**2) + self.h0) * np.sin(self.lat0)
        xyz0 = np.array([x0, y0, z0])
        xy0 = np.array([x0, y0])
        r0 = np.sqrt(np.dot(xyz0, xyz0))
        p0 = np.sqrt(np.dot(xy0, xy0))

        self.xyz0 = xyz0 # to be remembered between calls
        # unit vectors to be remembered between calls, hence self.xxx
        self.east0 = np.array([-y0, x0, 0]) / p0
        # zenith and north directions wrt local ellipsoid
        self.zenith0 = np.array([np.cos(self.lat0) * np.cos(self.lon0),
                                 np.cos(self.lat0) * np.sin(self.lon0),
                                 np.sin(self.lat0)])
        self.north0 = np.cross(self.zenith0, self.east0)


    def xyz2dec_ha(self, vec):
        """Declination and hour angle in target direction used to describe radar beam
        direction at JRO, corresponding to latitude and relative longitude of the
        beam-spot on the celestial sphere, corresponds to rr->\infty, in which case:
        """
        vec = vec / np.sqrt(np.dot(vec, vec))
        p = np.sqrt(vec[0]**2 + vec[1] ** 2)
        dec = np.arctan2(vec[2], p) / deg                                  # in degrees
        ha = (np.arctan2(vec[1], vec[0]) - self.lon0) * (24 / (2.*np.pi)) * 60  # in minutes

        return dec,ha

    def dec_ha2el_az(self, dec, ha):
        """Returns elevation and azimuth angles of a radar beam with respect to
        local tangent plane.

        the beam is specified by:
                declination dec (deg)
                hour angle  ha  (min)
        with respect to radar location at longitude lon0 and height h0
        above reference ellipsiod at geodetic latitude lat0
        """

        lat = dec * deg                                # on celestial sphere
        lon = 2 * np.pi * (ha / (24 * 60))
        lon = lon + lon0                    # on celestial sphere
        vec = np.array([np.cos(lat) * np.cos(lon),
                        np.cos(lat) * np.sin(lon),
                        np.sin(lat)])
        hor = vec - np.dot(vec, self.zenith0) * self.zenith0
        hor = hor / np.sqrt(np.dot(hor, hor))
        el = np.arccos(np.dot(hor, vec)) / deg
        north = np.dot(hor, self.north0)
        east = np.dot(hor, self.east0)
        az = np.arctan2(east, north) / deg

        return el,az

    def aspect_angle(self, year, xyz):
        """Returns the magnetic aspect angle (rad) of a target with
        geocentric vector xyz defined in geocentric coordinates
        """

        r = np.sqrt(np.dot(xyz, xyz))
        p = np.sqrt(xyz[0] ** 2 + xyz[1] ** 2)
        lat = np.arctan2(xyz[2], p)
        lon = np.arctan2(xyz[1], xyz[0])
        radial = xyz / r;                        # directions from target
        east = np.array([-xyz[1], xyz[0], 0]) / p
        north = -np.cross(east, radial)
        rr = xyz - self.xyz0
        u_rr = rr / np.sqrt(np.dot(rr, rr))      # unit vector from radar to target

        [bX, bY, bZ, bB] = self.igrf.igrf_B(year, r-a_igrf, lon / deg, lat / deg)
        bfield = np.array([bX, bY, bZ])
        B = bX * north + bY * east - bZ * radial
        u_B = B / np.sqrt(np.dot(B, B))
        aspect = np.arccos(np.dot(u_B, u_rr))
        return r, lat, lon, aspect


    def aspect_elaz(self, year, rr, el, az):
        """Returns magnetic aspect angle and geocentric coordinates of a target tracked by jro at
        range       rr (km)
        elevation   el (rad above local tangent plane to ellipsoid)
        azimuth     az (rad east of local north)
        """

        tx = np.cos(el) * np.sin(az)                    # direction cosines wrt east and north
        ty = np.cos(el) * np.cos(az)
        tz = np.sin(el)
        xyz = self.xyz0 + rr * (tx * self.east0 + ty * self.north0 + tz * self.zenith0)
        #geocentric coordinates of target

        [r, lat, lon, aspect] = self.aspect_angle(year, xyz)
        [dec,ha] = self.xyz2dec_ha(xyz - self.xyz0)
        return r, lon, lat, dec, ha, aspect

    def cosBs(self, year, rr, el, az):
        """Decomposes the radial unit vector to the target to direction cosines
        of magnetic North, East, and Up
        """

        tx = np.cos(el) * np.sin(az)                              # direction cosines wrt east and north
        ty = np.cos(el) * np.cos(az)
        tz = np.sin(el)
        xyz = self.xyz0 + rr * (tx * self.east0 + ty * self.north0 + tz * self.zenith0)     # target vector
        r = np.sqrt(np.dot(xyz, xyz))
        lat, lon, h = self.xyz2llh(xyz[0], xyz[1], xyz[2])         # target lat, lon, height

        radial = xyz / r; # unit vector to target
        p = np.sqrt(xyz[0] ** 2 + xyz[1] ** 2)
        east = np.array([-xyz[1], xyz[0], 0]) / p # unit vector to east from target
        north = -np.cross(east, radial) # unit vector to north from target
        rr_ = xyz - self.xyz0 # vector from radar to target
        rr_u = rr_ / np.sqrt(np.dot(rr_, rr_)) # unit vector from radar to target

        [bX, bY, bZ, bB] = self.igrf.igrf_B(year, r - a_igrf, lon / deg, lat / deg)
        bfield = np.array([bX, bY, bZ])
        B = bX * north + bY * east - bZ * radial # magnetic field vector B
        bn = B / np.sqrt(np.dot(B, B))
        # "magnetic north" unit vector since B points by definition in "magnetic north" direction

        be = np.cross(bn, radial)
        be = be / np.sqrt(np.dot(be, be)) # magnetic east unit vector
        bu = np.cross(be, bn) # magnetic up unit vector

        cosBn = np.dot(bn, rr_u) # magnetic north direction-cosine of rr_u
        aspect_angle = np.arccos(cosBn)
        cosBe = np.dot(be, rr_u) # magnetic east direction-cosine of rr_u
        cosBu = np.dot(bu, rr_u) # magnetic up direction-cosine of rr_u

        """
        uLOS=cosBe*U(h)+cosBn*V(h)+cosBu*W(h) ... LOS wind model in terms of wind components to calculate and direction cosines
        """

        return r, lat, lon, h, xyz, B, aspect, cosBn, cosBe, cosBu
