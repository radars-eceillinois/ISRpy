"""
#   jrobeam.py
#
#    module for calculating geometry parameters and magnetic aspect
#    angle of radar targets monitored by jro
#
#    use aspect_elaz or aspect_txty to calculate aspect angles of targets
#    specified by (el,az) or (tx,ty) angles
#
#  Created by Erhan Kudeki on 11/29/08.
#  Copyright (c) 2008 ECE, UIUC. All rights reserved.
#
"""

from .beam import *
import numpy as np

# ------------ jro radar specifications -------------------------
jrospecs = RadarSpecs(
    lat0 = -11.947917 * deg,   # geodetic, the usual map or GPS latitude
    lon0 = -76.872306 * deg,   # east of Greenwich
    h0 = 0.463,                # local height above reference ellipsoid
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = jrospecs.lat0    # geodetic, the usual map or GPS latitude
lon0 = jrospecs.lon0    # east of Greenwich
h0   = jrospecs.h0      # local height above reference ellipsoid


n0 = jrospecs.n0
x0 = jrospecs.x0
y0 = jrospecs.y0
z0 = jrospecs.z0
xyz0 = jrospecs.xyz0
uo = jrospecs.uo        # on axis

# unit vectors from jro
east0 = jrospecs.east0
zenith0 = jrospecs.zenith0
north0 = jrospecs.north0

# radar methods from beam.py
xyz2dec_ha = jrospecs.xyz2dec_ha
dec_ha2el_az = jrospecs.dec_ha2el_az
aspect_angle = jrospecs.aspect_angle
aspect_elaz = jrospecs.aspect_elaz
cosBs = jrospecs.cosBs

# orthonormal basis vectors including the jro on-axis direction
dec = -12.88 * deg,        # antenna declination
ha = -(4+37./60.) * deg    # hour angle on-axis direction at JRO

uo = np.array([np.cos(dec) * np.cos(ha/4. + lon0),
               np.cos(dec) * np.sin(ha/4. + lon0),
               np.sin(dec)])    # on axis
ux = np.cross(zenith0, uo)
ux = ux / np.sqrt(np.dot(ux, ux))  # along the building to the right
uy = np.cross(uo, ux)              # away from the building into the valley

def aspect_txty(year, rr, tx, ty):
    """Returns magnetic aspect angle and geocentric coordinates of a target tracked by jro at
    range rr (km)
    tx along jro building
    ty into the building
    """

    tz = np.sqrt(1 - tx ** 2. - ty ** 2.)
    xyz = xyz0 + rr * (tx * ux + ty * uy + tz * uo)
    #geocentric coordinates of target

    [r, lat, lon, aspect] = self.aspect_angle(year, xyz)
    [dec, ha] = xyz2dec_ha(xyz - xyz0)

    return r,lon,lat,dec,ha,aspect
