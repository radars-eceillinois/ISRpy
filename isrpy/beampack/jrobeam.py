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
    dec = -12.88 * deg,        # antenna declination
    ha = -(4+37./60.) * deg    # hour angle on-axis direction at JRO
    )

# Make this variables accesible to the module for backwards compatibility

lat0 = jrospecs.lat0    # geodetic, the usual map or GPS latitude
lon0 = jrospecs.lon0    # east of Greenwich
h0   = jrospecs.h0      # local height above reference ellipsoid
dec = jrospecs.dec      # antenna declination
ha  = jrospecs.ha       # hour angle on-axis direction at JRO

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
aspect_angle = jrospecs.aspect_angle
aspect_txty = jrospecs.aspect_txty
aspect_elaz = jrospecs.aspect_elaz
cosBs = jropecs.cosBs

# orthonormal basis vectors including the jro on-axis direction

ux = np.cross(zenith0, uo)
ux = ux / np.sqrt(np.dot(ux, ux))  # along the building to the right
uy = np.cross(uo, ux)              # away from the building into the valley
